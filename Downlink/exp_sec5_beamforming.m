function exp_sec5_beamforming()
close all
%% ========  Pre SET Parameters  ======= %
% Read EVM file
dir          = pwd;
N            = 3; % Number of antennas
fc           = 915e6;
c_light      = 3e8;
lambda       = c_light/fc;
dx           = lambda/2;
%
%% ============= STEP 1 ================ %
BW = 125e3; Fs = BW;
data = zeros(N,2.12e5);
for select_mode = 1:3
if     select_mode == 1
    str = '_123456789#';
    weight_factor = 0.3;
    angle_deg      = -20;

elseif select_mode == 2
    str = '_ABCDEFGHJKL';
    weight_factor = 0.4 ;
    angle_deg      = 30;
elseif select_mode == 3
    str = '_HelloWorld';
    angle_deg      = 70;
    weight_factor = 0.3;
elseif select_mode == 4
    str = '_listening4x4';
    angle_deg      = 90;
end

[raw_data,EXPORT_ROOT] = GetDir_Read_lora_Data(str);
% figure;plot(abs(raw_data_1to9),'k');hold on; plot(abs(raw_a2l),'r');
% FigTFSpectrum([raw_a2l(1:length(raw_data_1to9))+raw_data_1to9],Fs);
% [raw_data_1to9,EXPORT_ROOT] = GetDir_Read_lora_Data('_123456789#');
% [raw_a2l,EXPORT_ROOT] = GetDir_Read_lora_Data('_ABCDEFGHJKL');
% [raw_helloworld,EXPORT_ROOT] = GetDir_Read_lora_Data('_HelloWorld');
% index    = find(raw_data_1to9~=0);
% raw_data = [zeros(index(1)+100,1);raw_data]; 
% figure;plot(abs(raw_data_1to9),'k');hold on; plot(abs(raw_a2l),'r');hold on;
%  plot(abs(raw_helloworld),'b');hold on;
% hold on; plot(abs(raw_data),'b');
theta = linspace(-pi/2, pi/2, 360); % Angles from -90 to 90 degrees in radians
% FigTFSpectrum([raw_a2l(1:length(raw_data_1to9))+raw_data_1to9],Fs);
% Steering angle (incoming signal direction)
steering_angle = angle_deg * pi / 180; % 30 degrees converted to radians

% Steering vector
steering_vector = exp(1j * (0:(N-1))' * (2*pi*dx/lambda) * sin(theta - steering_angle));

% Apply Hamming window to the steering vector
window = hanning(N);
windowed_steering_vector = window .* steering_vector;

% Array factor (Beamforming)
AF = sum(windowed_steering_vector, 1);

% Normalize the array factor
AF_normalized = AF / max(abs(AF));

% Convert to dB scale
AF_dB = 20 * log10(abs(AF_normalized));
% Plot the antenna pattern in polar coordinates
figure;
polarplot(theta, AF_dB);
title('Four-Antenna Beamforming Pattern at 915 MHz with Hamming Window');
rticks([-40 -30 -20 -10 0]); % Define radius ticks in dB
rlim([-40 0]); % Define radius limits
thetalim([-90 90]); % Define theta limits from -90 to 90 degrees
set(gca, 'ThetaZeroLocation', 'top'); % Set zero angle to top
set(gca, 'ThetaDir', 'clockwise'); % Set angle direction to clockwise
grid on;
%% ============= STEP 2 Beamforming ================ %
beamforming_vector = window .*exp(1j * (0:(N-1))' * (2*pi*dx/lambda) * sin(steering_angle));
% apply weights to the signal
for ii = 1 : N
    weight_sig(ii,:) = weight_factor* beamforming_vector(ii) * [raw_data];
end
  data =   weight_sig+data;
end
 strr = 'downlink';
if N ==4 
  EXPORT_ROOT = [EXPORT_ROOT , '4_antenna/'];
elseif N ==3 
     EXPORT_ROOT = [EXPORT_ROOT , '3_antenna/'];
end
for id = 1 : N
write_Combpkt_fileSource(data(id,:),EXPORT_ROOT,strr,id);
end
end



%=================================================================%
% %                sub function
%=================================================================%
function [raw_data,EXPORT_ROOT] = GetDir_Read_lora_Data(str);
FileHead = 'RewriteSource';
Filestr = [FileHead,str];
DataDIR  = '/RAWDATA/LoRa_CtrlPKT/Tx_packet/';
dir      = pwd;
FileDIR  = [dir,DataDIR];
INPORT_ROOT  = [FileDIR];
%  Data save path
FileName     = [INPORT_ROOT, Filestr '.dat'];
fi_1         = fopen(FileName, 'rb');
values       = fread(fi_1, Inf, 'float32');
raw_data     = values(1:2:end) + 1j*values(2:2:end);
fclose(fi_1);
disp(['Read the ', Filestr,' data']);
EXPORT_ROOT = [dir,DataDIR];
end
%-------------------------------------------------------------------------------------%
%-------------------------------------------------------%

function write_Combpkt_fileSource(sig_pkt,EXPORT_ROOT,str, tx_id)
% This function is to write new packet

DIR =  [EXPORT_ROOT];
file_sig = fopen([DIR, 'Beamforming_', str, 'Tx',num2str(tx_id),'.dat'], 'w');
sig_data = zeros(2*length(sig_pkt), 1);
sig_data(1:2:end) = real(sig_pkt);
sig_data(2:2:end) = imag(sig_pkt);
fwrite(file_sig, sig_data, 'float32');
fclose(file_sig);
disp(['Write the ', str,num2str(tx_id),' data...']);

end
%-----------------------------------------------------%
%-----------------------------------------------------%

function [fig] = FigTFSpectrum(signal,Fs,varargin);
%   Function - SPECTRUM

if (nargin < 2)
    error(message('Error'));
end
if     (nargin == 2)
    win     = 32;
    nfft    = 256;
    LOG_PSD = 1;
    str     = 'TF Spectrum';
elseif (nargin == 3)
    str     = varargin{1};
    win     = 32;
    nfft    = 256;
    LOG_PSD = 1;
elseif (nargin == 4)
    str     = varargin{1};
    win     = varargin{2};
    nfft    = 256;
    LOG_PSD = 1;
elseif (nargin == 5)
    str     = varargin{1};
    win     = varargin{2};
    nfft    = varargin{3};
    LOG_PSD = 1;
elseif (nargin == 6)
    str     = varargin{1};
    win     = varargin{2};
    nfft    = varargin{3};
    LOG_PSD = varargin{4};
end

figure;
overlap = floor(win*7/8);
[S,F,T,P] = spectrogram(signal, win, overlap, nfft, Fs, 'centered');
norm_freq = F/1000;     % unit: kHz
win_id = 1:1:numel(T);  % window ID
if (LOG_PSD == 1)
    psd = 10*log10(abs(P));
else
    psd = abs(P);
end
fig = surf(win_id, norm_freq, psd);
axis tight;
shading interp;
view(0, 90);
ylim([-1 1]*max(norm_freq));
title(str)
xlabel('Win ID');
ylabel('Frequency (kHz)');
grid off;
%colorbar
set(gcf, 'unit', 'centimeters', 'position', [1/4, 8, 40, 5]);
left_margin = 0.05;
right_margin = 0.02;
bot_margin = 0.22;
top_margin = 0.15;
set(gca, 'position', [left_margin, bot_margin, ...
    1-left_margin-right_margin, 1-bot_margin-top_margin]);

end
%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%





