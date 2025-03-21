function sec4_beamforming(N_antenna, angle_deg, select_mode)
close all
%% ========  Pre SET Parameters  ======= %
% Read EVM file
dir          = pwd;
fc           = 915e6;
c_light      = 3e8;
lambda       = c_light/fc;
dx           = lambda/2;

%% ============= STEP 1 ================ %
BW = 125e3; Fs = BW;

if     select_mode == 1
    str = '123456789#';
elseif select_mode == 2
    str = 'ABCDEFGHJKL';
elseif select_mode == 3
    str = 'HelloWorld';
elseif select_mode == 4
    str = 'listening4x4';
end

[raw_data,EXPORT_ROOT] = GetDir_Read_lora_Data(str);
% figure;plot(abs(raw_data_1to9),'k');hold on; plot(abs(raw_a2l),'r');
% FigTFSpectrum([raw_a2l(1:length(raw_data_1to9))+raw_data_1to9],Fs);
% [raw_data_1to9,EXPORT_ROOT] = GetDir_Read_lora_Data('1to9');
% [raw_a2l,EXPORT_ROOT] = GetDir_Read_lora_Data('A2L');
% %  figure;plot(abs(raw_data_1to9),'k');hold on; plot(abs(raw_a2l),'r');
theta = linspace(-pi/2, pi/2, 360); % Angles from -90 to 90 degrees in radians

% Steering angle (incoming signal direction)
steering_angle = angle_deg * pi / 180; % 30 degrees converted to radians
% Steering vector
steering_vector = exp(1j * (0:(N_antenna-1))' * (2*pi*dx/lambda) * sin(theta - steering_angle));

% Apply Hamming window to the steering vector
window = hanning(N_antenna);
windowed_steering_vector = window .* steering_vector;

% Array factor (Beamforming)
AF = sum(windowed_steering_vector, 1);

% Normalize the array factor
AF_normalized = AF / max(abs(AF));

% Convert to dB scale
AF_dB = 20 * log10(abs(AF_normalized));
% Plot the antenna pattern in polar coordinates
% figure;
% polarplot(theta, AF_dB);
% title('Four-Antenna Beamforming Pattern at 915 MHz with Hamming Window');
% rticks([-40 -30 -20 -10 0]); % Define radius ticks in dB
% rlim([-40 0]); % Define radius limits
% thetalim([-90 90]); % Define theta limits from -90 to 90 degrees
% set(gca, 'ThetaZeroLocation', 'top'); % Set zero angle to top
% set(gca, 'ThetaDir', 'clockwise'); % Set angle direction to clockwise
% grid on;
%% ============= STEP 2 Beamforming ================ %
beamforming_vector = window .*exp(1j * (0:(N_antenna-1))' * (2*pi*dx/lambda) * sin(steering_angle));
% apply weights to the signal
for ii = 1 : N_antenna
    weight_sig(ii,:) = beamforming_vector(ii) * raw_data;
end

for id = 1 : N_antenna
write_pkt_fileSource(weight_sig(id,:),EXPORT_ROOT,str,id);
end
end



%=================================================================%
% %                sub function
%=================================================================%
function [raw_data,EXPORT_ROOT] = GetDir_Read_lora_Data(str);
FileHead = 'RewriteSource_';
Filestr = [FileHead,str];
DataDIR  = '/RAWDATA/LoRa_CtrlPKT/';
dir      = pwd;
FileDIR  = [dir,DataDIR];
INPORT_ROOT  = [FileDIR,'input/'];
%  Data save path
FileName     = [INPORT_ROOT, Filestr '.dat'];
fi_1         = fopen(FileName, 'rb');
values       = fread(fi_1, Inf, 'float32');
raw_data     = values(1:2:end) + 1j*values(2:2:end);
fclose(fi_1);
disp(['Read the ', Filestr,' data']);
EXPORT_ROOT = [dir,'/RAWDATA/USRP/', 'Tx_packet/'];
end
%-------------------------------------------------------------------------------------%
%-------------------------------------------------------%

function write_pkt_fileSource(sig_pkt,EXPORT_ROOT,str, tx_id)
% This function is to write new packet

DIR =  [EXPORT_ROOT];
file_sig = fopen([DIR, 'Beamforming_', str, 'Tx',num2str(tx_id),'.dat'], 'w');
sig_data = zeros(2*length(sig_pkt), 1);
sig_data(1:2:end) = real(sig_pkt);
sig_data(2:2:end) = imag(sig_pkt);
fwrite(file_sig, sig_data, 'float32');
fclose(file_sig);
%disp(['Write the ', str,num2str(tx_id),' data...']);

end
%-----------------------------------------------------%

%--------------------------------------------------------------------------------------%





