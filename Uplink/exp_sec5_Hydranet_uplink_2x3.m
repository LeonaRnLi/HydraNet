function exp_sec5_Hydranet_uplink_3x2()
% ---------------------------------------------------------------------------- %
%%  This function is to do the pre-process module
% STEP 1: extract the packets
% ---------------------------------------------------------------------------- %
%
close all

%% ========  Pre SET Parameters  ======= %
% Read EVM file
dir          = pwd;
Dir_all      = [dir, '/RAWDATA/LoRa_Tx_packet/'];
str          = 'PurePacket/';
realFileDIR  = [Dir_all, str, 'input/'];
[SF,BW,Fs]   = Parameter_xml_Document_Read();
FileHead = {'123456789#','ABCDEFGHJKL','HelloWorld','listening4x4'};
Filename1= ['normal_demod' 'Payload_NODE_', FileHead{1}];
Filename2= ['normal_demod' 'Payload_NODE_', FileHead{2}];
Filename3= ['normal_demod' 'Payload_NODE_', FileHead{3}];

file_info1   = [realFileDIR, Filename1,'.csv'];
file_info2   = [realFileDIR, Filename2,'.csv'];
file_info3   = [realFileDIR, Filename3,'.csv'];


packet1      = csvread(file_info1);
packet2      = csvread(file_info2);
packet3      = csvread(file_info3);
%% ============= STEP 1 ================ %
%  Extract the two channel packets
%  we will not estimate the first packet.
%  we can choose the second or third.
[raw_data1,EXPORT_ROOT] = GetDir_Read_lora_Data(1);
[raw_data2]             = GetDir_Read_lora_Data(2);
% [raw_data3]             = GetDir_Read_lora_Data(3);
% [raw_data4]             = GetDir_Read_lora_Data(4);

RX_ID1 = 1;RX_ID2 = 2;RX_ID3 = 3;RX_ID4 = 4;
chirplen = 2^SF ;

threshold   = 5e-2; preMaxRange   = 2^SF * 200; pkt_len =  chirplen*3;
FRMlen      = Fs * 10;CoarseStep  = chirplen* 10;
display_on  = 0;

% extract_signalframe_location(raw_data1,EXPORT_ROOT,threshold,display_on, pkt_len, preMaxRange,CoarseStep,FRMlen,1);
% extract_signalframe_location(raw_data2,EXPORT_ROOT,threshold,display_on, pkt_len, preMaxRange,CoarseStep,FRMlen,2);
% % % % % extract_signalframe_location(raw_data3,EXPORT_ROOT,threshold,display_on, pkt_len, preMaxRange,CoarseStep,FRMlen,3);
% % % % % extract_signalframe_location(raw_data4,EXPORT_ROOT,threshold,display_on, pkt_len, preMaxRange,CoarseStep,FRMlen,4);

%  Read the cluster
Name          = 'frameLocs_cluster_Rx';
cluster_info1 = csvread([EXPORT_ROOT [Name, num2str(RX_ID1), '.csv']]);
cluster_info2 = csvread([EXPORT_ROOT [Name, num2str(RX_ID2), '.csv']]);
% cluster_info3 = csvread([EXPORT_ROOT [Name, num2str(RX_ID3), '.csv']]);
% cluster_info4 = csvread([EXPORT_ROOT [Name, num2str(RX_ID4), '.csv']]);
n_frame = size(cluster_info1,1)-1;
startnum1   = cluster_info1(1,2);
startnum2   = cluster_info1(1,2);
% startnum3   = cluster_info4(1,2);
% startnum4   = cluster_info4(1,2);
N_antenna   = 2;
weight_vary = [0.25,0.5,1:2:6];

%% ============= STEP 2 ================ %
%  Extract a packet
for pkt_id = 1:n_frame
    for ww = 1:length(weight_vary)
        vary = 2^SF*50;
        
        relativeID     = pkt_id +1;
        cluster_st1     = startnum1 + cluster_info1(relativeID,1);
        cluster_ed1     = cluster_st1 + vary;
        IQ_sample1      = raw_data1(cluster_st1 : cluster_ed1);
        cluster_st2     = startnum2 + cluster_info1(relativeID,1);
        cluster_ed2     = cluster_st1 + vary;
        IQ_sample2      = raw_data2(cluster_st2 : cluster_ed2);
       
        Rx_data = [IQ_sample1,weight_vary(ww)*IQ_sample2];
    
        deg_vary  = -90:5:90;
        
        for angle = 1:length(deg_vary)
            angle_in = deg_vary(angle);
            
            [beamforming_vector2] = generate_beam_weight(N_antenna,angle_in);
            Rx_post_decode2 = Rx_data*beamforming_vector2;
            [num_error21(ww,angle),num_error22(ww,angle),num_error23(ww,angle)] ...
                = compare_calculate_ser(Rx_post_decode2,packet1,packet2,packet3);
        end
    end

    [n_1to9,idx_1to9] = min(min([num_error21]));
    %     angle1(pkt_id) = deg_vary(idx_1to9);
    [n_a2l, idx_a2l]  = min(min([num_error22]));
    %     angle2(pkt_id) = deg_vary(idx_a2l);
    
    [n_hw, idx_hw]    = min(min([num_error23]));
    %     angle3(pkt_id) = deg_vary(idx_hw);
    error_num(pkt_id,:) = [n_1to9,n_a2l,n_hw];
    disp(['Error num: ',num2str(error_num(pkt_id,1)),' ',num2str(error_num(pkt_id,2)),' ',num2str(error_num(pkt_id,3))]);
    filename_str = 'experiment_Hydranet_3x2';
    disp('save')
    append_to_excel(filename_str, error_num(pkt_id,:))
end
% filename_str = 'experiment_3x4';
% append_to_excel(filename_str, error_num(pkt_id,:))

end





%-------------------------------------------------------------------------------------%
%%                            SUB- FUNCTION
%-------------------------------------------------------------------------------------%
function [num_error11,num_error12,num_error13] = compare_calculate_ser(Rx_post_decode1,packet1,packet2,packet3)
[message1] = LoRa_Demodulate_normal(Rx_post_decode1,2);
if message1 ~= 0
    [n_ber1_1to9] = calculate_ber(message1,packet1);
    [n_ber1_A2L]  = calculate_ber(message1,packet2);
    [n_ber1_HW]   = calculate_ber(message1,packet3);
    num_error11 = min(n_ber1_1to9);
    num_error12 = min(n_ber1_A2L);
    num_error13 = min(n_ber1_HW);
else
%     disp('Not Detect!');
    num_error11 = 55;
    num_error12 = 55;
    num_error13 = 55;
end
end
%----------------------------------------------------------------------%
function append_to_excel(filename, data)
% 检查文件是否存在
path = [pwd, '/SystemPara/', 'RAWDATA/',filename,'.csv'];
% 检查文件是否存在
if exist(path, 'file') == 2
    % 如果文件存在，读取现有数据
    existing_data = dlmread(path, ',');
else
    % 如果文件不存在，设置 existing_data 为空数组
    existing_data = [];
end

% 将现有数据和新数据合并
combined_data = [existing_data; data];

% 写入新数据
dlmwrite(path, combined_data,'precision', 9);

end

%-------------------------------------------------------------------------------------%

function [beamforming_vector] = generate_beam_weight(N_antenna,angle_deg)
theta = linspace(-pi, pi, 360); % Angles from -90 to 90 degrees in radians
fc           = 915e6;
c_light      = 3e8;
lambda       = c_light/fc;
dx           = lambda/2;

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
% title(['Beamforming Degree: ',num2str(angle_deg)]);
% rticks([-40 -30 -20 -10 0]); % Define radius ticks in dB
% rlim([-40 0]); % Define radius limits
% thetalim([-90 90]); % Define theta limits from -90 to 90 degrees
% set(gca, 'ThetaZeroLocation', 'top'); % Set zero angle to top
% set(gca, 'ThetaDir', 'clockwise'); % Set angle direction to clockwise
% grid on;
%% ============= STEP 2 Beamforming ================ %
beamforming_vector = window .*exp(1j * (0:(N_antenna-1))' * (2*pi*dx/lambda) * sin(steering_angle));



end
%-------------------------------------------------------------------------------------%
function [raw_data,EXPORT_ROOT] = GetDir_Read_lora_Data(flag);
FileHead = 'USRP_PostCoding_LoRa_RX';
Filestr = [FileHead,num2str(flag)];
DataDIR  = '/RAWDATA/USRP/Hydranet_LoRa/';
dir      = pwd;
FileDIR  = [dir,DataDIR];
INPORT_ROOT  = [FileDIR,'input_3x3/'];
EXPORT_ROOT  = [FileDIR, 'output_3x3/Cluster/'];
%  Data save path
FileName     = [INPORT_ROOT, Filestr '.dat'];
fi_1         = fopen(FileName, 'rb');
values       = fread(fi_1, Inf, 'float32');
raw_data     = values(1:2:end) + 1j*values(2:2:end);
fclose(fi_1);
disp(['Read the ', Filestr,' data']);
end
%-------------------------------------------------------------------------------------%
function [n_ber,ber] = calculate_ber(message,int_pkt);
num_pkt1   = length(message)-length(int_pkt)+1;
if num_pkt1<=0
    num_pkt1 = 1;
    int_pkt = int_pkt(1:length(message));
end
for ii = 1 : num_pkt1
    error1 = abs(message(ii+[0:length(int_pkt)-1])-int_pkt);
    n_ber(ii) = length(find(error1>3));
end
ber = n_ber/length(int_pkt);

end

%-------------------------------------------------------------------------------------%
%%=====================================================%%
%%          Function - Main- demod
%%=====================================================%%
function [symbols_message,symbols_Demod] =...
    LoRa_Demodulate_normal(signal,Coherece)
% This function is aim to demodulate a received signal
%  in the normal way.
%  input
%        signal       IQ LoRa signal containing (preamble + sync header + payload)
%        Coherence    type of demodulation (coherent or non-coherent)
%                     (1) coherent or (2) non-coherent FSK Detection
%        Fs           sampling fre1.quency
%        dF           carrier frequency shift (Fc-fc) % Fc means Sampling
%                     central freq, fc means carrier freq
%
%  output
%        symbols_message          message symbols vector (encoded)
%        symbols_Demod            LoRa symbols vector
%        n_preamble               Number of symbols in preamble
% Ruonan - Jun 23 2023
[SF, BW,Fs,~,FFT_factor] = Parameter_xml_Document_Read();
if Fs == BW
    SNR = Inf;
    signal_demod        = awgn(signal,SNR,'measured') ;
else
    df= 0; SNR= Inf ;
    signal_freq_demod   = signal.*exp(1j.*2.*pi.*df./Fs.*(0:length(signal)-1))' ;
    LP_filter           = fir1(40,BW/(Fs/2),'low',chebwin(41,30));
    signal_filter = filter(LP_filter,1,signal_freq_demod);
    signal_demod        = awgn(resample(signal_filter,BW,Fs),SNR,'measured') ;
end
%  Return if SF is not in the range
if SF > 12 || SF < 7
    return
end
chirplen  = 2^SF *FFT_factor;
signal    = signal_demod;
[MessageStartInd,symbol_offset, sfd_off] = Preamble_detect_Corr(signal);
if MessageStartInd == 0
    symbols_message = 0;
    return
end
Nmessage        = floor(length(signal)/chirplen - (MessageStartInd/chirplen)) ;
MessageEndInd   = Nmessage.*chirplen + MessageStartInd ;
MessageSignal   = signal(MessageStartInd+1:MessageEndInd).*...
    loramod((sfd_off)*ones(1,Nmessage),SF,BW,BW,-1) ;
symbols_Demod   = FSKDetection(MessageSignal,SF,Coherece) ;
symbols_message = mod(symbols_Demod - symbol_offset,chirplen);
end

%-----------------------------------------------------%
%   Function - multi dechirp signal generate
%-----------------------------------------------------%
function [y] = loramod(x,SF,BW,fs,varargin)
% loramod LoRa modulates a symbol vector specified by x
%
%   in:  x          1xN symbol vector wher N=1-Inf
%                   with values {0,1,2,...,2^(SF)-1}
%        BW         signal bandwidth of LoRa transmisson
%        SF         spreading factor
%        Fs         sampling frequency
%        varargin{1} set polarity of chirp
%
%  out:  y          LoRa IQ waveform
if (nargin < 4)
    error(message('comm:pskmod:numarg1'));
end

if (nargin > 5)
    error(message('comm:pskmod:numarg2'));
end

% Check that x is a positivFs = BW;e integer
if (~isreal(x) || any(any(ceil(x) ~= x)) || ~isnumeric(x))
    error(message('comm:pskmod:xreal1'));
end
M       = 2^SF ;
% Check that M is a positive integer
if (~isreal(M) || ~isscalar(M) || M<=0 || (ceil(M)~=M) || ~isnumeric(M))
    error(message('comm:pskmod:Mreal'));
end

% Check that x is within range
if ((min(min(x)) < 0) || (max(max(x)) > (M-1)))
    error(message('comm:pskmod:xreal2'));
end

% Polarity of Chirp
if nargin == 4
    Inv = 1 ;
elseif nargin == 5
    Inv = varargin{1} ;
end
% Symbol Constants
Ts      = 2^SF/BW ;
Ns      = fs.*M/BW ;

gamma   = x/Ts ;    %f0
beta    = BW/Ts ;   %Kr

time    = (0:Ns-1)'.*1/fs ;
freq    = mod(gamma + beta.*time,BW) - BW/2 ;
%  figure;plot(Inv.*freq(:,2)*2/BW)
Theta   = cumtrapz(time,Inv.*freq) ;
y       = reshape(exp(1j.*2.*pi.*Theta),numel(Theta),1) ;
end
%-----------------------------------------------------%
%   Function - PREAMBLE & SFD DETECT
%-----------------------------------------------------%
function [SFD_edge,prmb_offset,offset] =...
    Preamble_detect_Corr(signal,varagin)
% This main function is aim to
% detect preamble edge using correlation
[SF, BW] = Parameter_xml_Document_Read();
Fs = BW;
if     (nargin == 1)
    Coherece   = 2;
    pream_len  = 8;
    TOLERANCE  = 2;
elseif (nargin == 2)
    pream_len  = 8;
    TOLERANCE  = 2;
elseif (nargin == 3)
    TOLERANCE  = 2;
elseif (nargin < 1)
    return;
    disp('Error: Not input any signal!');
end

chirplen       = Fs/BW * 2^SF;

if  (length(signal) < pream_len * chirplen)
    return;
    disp('Error: signal is too short!');
end

Nsymbols        = floor(length(signal)/chirplen) ;
DChirpsDemod    = loramod(zeros(1,Nsymbols),SF,BW,BW,-1);
signal_demod    = signal(1:length(DChirpsDemod));
SniffSignal     = signal_demod.*DChirpsDemod ;
PreamInd        = FSKDetection(SniffSignal,SF,Coherece);
found = 0;
ii = 0;

while found == 0
    if (ii +  pream_len  >= length(PreamInd))
%         disp('Error: can not find preamble!')
        SFD_edge = 0;
        prmb_offset = 0;
        offset = 0;
        return;
        
    end
    pream_check  = PreamInd(ii + (1 : pream_len));
    idx          = pream_check(1);
    [found]      = Check_INDEXTolerance(idx,pream_check,TOLERANCE);
    ii           = ii + 1;
end

% disp('Detect preamble, SUCCESS!');
prmb_offset      = mode(pream_check);
%-----------------------------------------------------%
UChirpsDemod     = loramod(zeros(1,Nsymbols),SF,BW,BW);
signal_demod     = signal(1:length(UChirpsDemod));
SniffSignal      = signal_demod.*UChirpsDemod ;
fft_signal       = fft(reshape(SniffSignal,chirplen,...
    length(SniffSignal)/chirplen));
[~,SFDInd]       = sort(max(abs(fft_signal)));
SFD_rough        = sort(SFDInd(end-1:end)) ;
SyncRough        = SFD_rough(end) + 1  - 5;
%
st         = SyncRough * chirplen;
winlen     = 2^4;
SFD_search = [];
found      = 0;ii = 0;
while (found==0)
    %     st
    if (st +2 * chirplen >= length(signal))
%         disp('Error: can not find sfd!')
        SFD_edge = 0;
        prmb_offset = 0;
        offset = 0;
        return;
    end
    datain = signal(st + (1 : 2 * chirplen));
    for of = 0 : 1 : ( winlen - 1)
        offset     = fix(mod((of * chirplen / winlen) , chirplen));
        temp       = datain(offset + (1:chirplen));
        off_chirp  = loramod(offset,SF,BW,BW);
        [max_idx]  = FSKDetection(temp .* off_chirp,SF,2);
        SFD_search = [max_idx,SFD_search];
        idx        = SFD_search(1);
        if length(SFD_search) > 2 * winlen
            SFD_search(end)   = [];
            found = Check_INDEXTolerance(idx,SFD_search,TOLERANCE*1);
        end
        if (found)
            SFD_offset = fix(mod((offset + (chirplen / 4)) , chirplen));
            break;
        end
    end
    if ii > winlen
        
%         disp(['can not find SFD']) ;
        SFD_offset = 0;
        break;
    end
    st = st + chirplen;
end
SFD_edge = st + SFD_offset;
% disp('Detect SFD, SUCCESS!');

end

%-----------------------------------------------------%
%   Function - SYMBOL INITIAL FREQ
%-----------------------------------------------------%

function [symbols] = FSKDetection(signal,SF,detection)
% LoRa_Tx demodulates a Lora de-chirped signal using
% the coherence specified by the detection variable
%
%   in:  message      payload message
%        SF           spreading factor
%        detection    1= coherent detection, 2= non-coherent detection
%
%  out:  symbols      FSK demodulated symbol vector
[~, ~,~,~,fft_factor] = Parameter_xml_Document_Read();
chirp_len  = 2^SF;
num_fft    =  fft_factor * chirp_len;

if (detection == 1) % coherent detection
    t  = 0 : 1/(2^SF) : 0.999 ; % time vector
    for Ctr   = 1 : 2^SF
        rtemp = conv(signal,exp(-1j.*2.*pi.*(2^SF - Ctr + 1).*t)) ;
        % convolution w/ideal fsk signal
        r(Ctr,:) = real(rtemp(2^SF+1:2^SF:end)) ; % save resultant array
    end
    [~,idx]  = max(r) ; % take max
    symbols  = idx - 1 ; % store symbol vector
elseif (detection  == 2) % non-coherent detection
    signal_reshape = reshape(signal,2^SF,length(signal)/(2^SF));
    fft_result = fft(signal_reshape,num_fft);
    %     figure;imagesc(abs(signal_reshape))
    %     figure(4);plot(abs(fft_result(:,1)));hold on
    [~,idx] = max((fft_result)) ; % take max of fft window
    symbols = idx - 1 ; % store symbol array
end
end
%-----------------------------------------------------%
%   Function - Check index TOLERANCE
%-----------------------------------------------------%
function [found] = Check_INDEXTolerance...
    (idx,idx_history,TOLERANCE)
% INPUT:
%      idx:          standard idx
%      history:      compared idx history
%      rTOLERANCE:   limitation length
% OUTPUT:
%      found:        state flag (bool)
looplen = length(idx_history);
found   = 1;
if idx>= 254
    idx = abs(idx-256);
end
for ii  = 1 : 1 : (looplen - 1)
    if idx_history(1 + ii)>= 254
        idx_history(1 + ii) = abs( idx_history(1 + ii)-256);
    end
    
    if (abs(idx - idx_history(1 + ii)) > TOLERANCE)
        found = 0;
    end
end
end


%-----------------------------------------------------%
%   Function - SPECTRUM
%-----------------------------------------------------%

function [fig] = FigTFSpectrum(signal,Fs,varargin);

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


