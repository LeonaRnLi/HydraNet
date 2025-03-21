function exp_sec5_Hydranet_uplink_2x2()
% ---------------------------------------------------------------------------- %
%%  This function is to do the pre-process module
% STEP 1: extract the two channel packets
% STEP 2: do the packet detection, cfo estimate and fine time estimation
% STEP 3: estimate the channel
% STEP 4: save the estimated channel information.
% STEP 5: precode...and save.
%
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
% Filename3= ['normal_demod' 'Payload_NODE_', FileHead(3)];

file_info1   = [realFileDIR, Filename1,'.csv'];
file_info2   = [realFileDIR, Filename2,'.csv'];


packet1      = csvread(file_info1);
packet2      = csvread(file_info2);

%% ============= STEP 1 ================ %
%  Extract the two channel packets
%  we will not estimate the first packet.
%  we can choose the second or third.
[raw_data1,EXPORT_ROOT] = GetDir_Read_lora_Data(1);
[raw_data2]             = GetDir_Read_lora_Data(2);
[raw_data3]             = GetDir_Read_lora_Data(3);
% [raw_data4]             = GetDir_Read_lora_Data(4);

RX_ID1 = 1;RX_ID2 = 2;RX_ID3 = 3;
chirplen = 2^SF ;

threshold   = 5e-2; preMaxRange   = 2^SF * 200; pkt_len =  chirplen*3;
FRMlen      = Fs * 10;CoarseStep  = chirplen* 10;
display_on  = 1;
extract_signalframe_location(raw_data1,EXPORT_ROOT,threshold,display_on, pkt_len, preMaxRange,CoarseStep,FRMlen,1);
extract_signalframe_location(raw_data2,EXPORT_ROOT,threshold,display_on, pkt_len, preMaxRange,CoarseStep,FRMlen,2);
extract_signalframe_location(raw_data3,EXPORT_ROOT,threshold,display_on, pkt_len, preMaxRange,CoarseStep,FRMlen,3);
%  Read the cluster
Name          = 'frameLocs_cluster_Rx';
cluster_info1 = csvread([EXPORT_ROOT [Name, num2str(RX_ID1), '.csv']]);
cluster_info2 = csvread([EXPORT_ROOT [Name, num2str(RX_ID2), '.csv']]);
cluster_info3 = csvread([EXPORT_ROOT [Name, num2str(RX_ID3), '.csv']]);
% cluster_info4 = csvread([EXPORT_ROOT [Name, num2str(RX_ID4), '.csv']]);
n_frame = size(cluster_info1,1)-1;
startnum1   = cluster_info1(1,2);
startnum2   = cluster_info1(1,2);
startnum3   = cluster_info1(1,2);
% startnum3   = cluster_info4(1,2);
%% ============= STEP 2 ================ %
%  Extract a packet
for pkt_id = 1:n_frame
    
    vary = 2^SF*50;

    relativeID     = pkt_id +1;
    cluster_st1     = startnum1 + cluster_info1(relativeID,1);
    cluster_ed1     = cluster_st1 + vary;
    IQ_sample1      = raw_data1(cluster_st1 : cluster_ed1);
    cluster_st2     = startnum2 + cluster_info1(relativeID,1);
    cluster_ed2     = cluster_st1 + vary;
    IQ_sample2      = raw_data2(cluster_st2 : cluster_ed2);
    cluster_st3     = startnum3 + cluster_info1(relativeID,1);
    cluster_ed3     = cluster_st1 + vary;
    IQ_sample3      = raw_data3(cluster_st3 : cluster_ed3);
    % cluster_st4     = startnum4 + cluster_info4(relativeID,1);
    % cluster_ed4     = startnum4 + cluster_info4(relativeID,2);
    % IQ_sample4      = raw_data4(cluster_st4 : cluster_ed4);
    Rx_data = [IQ_sample1,IQ_sample2,IQ_sample3];
    % figure;plot(abs(IQ_sample1),'b');hold on ;plot(abs(IQ_sample2),'k');
    % title('Extract packet');xlabel('Sample ID');ylabel('Amplitude');legend('Rx1','RX2')
    
%     FigTFSpectrum(IQ_sample1,Fs,'RX1 Received Packet TF Spectrum');
%     FigTFSpectrum(IQ_sample2,Fs,'RX2 Received Packet TF Spectrum');
    % FigTFSpectrum(IQ_sample3,Fs,'RX3 Received Packet TF Spectrum');
    % FigTFSpectrum(IQ_sample3,Fs,'RX4 Received Packet TF Spectrum');
    angle_deg1 = -30;
    
    N_antenna= 3;
    [beamforming_vector1] = generate_beam_weight(N_antenna,angle_deg1);
%     [beamforming_vector2] = generate_beam_weight(N_antenna,angle_deg2);
    Rx_post_decode1 = Rx_data*beamforming_vector1;
%     Rx_post_decode2 = Rx_data*beamforming_vector2;
%     FigTFSpectrum(Rx_post_decode1,Fs,'RX1 Received Packet TF Spectrum');
%    deg_vary  = -90:10:90;
%    for angle = 1:length(deg_vary)
%        angle_in = deg_vary(angle);
angle_in = 35;
         [beamforming_vector2] = generate_beam_weight(N_antenna,angle_in);
         Rx_post_decode2 = Rx_data*beamforming_vector2;
%     FigTFSpectrum(Rx_post_decode2,Fs,['RX2 Received Packet angle ',num2str(angle_in)]);
%    end
%     
    
    
    
    [message1] = LoRa_Demodulate_normal(Rx_post_decode1,2);
    if message1 ~= 0
        [n_ber1_1to9] = calculate_ber(message1,packet1);
        [n_ber1_A2L] = calculate_ber(message1,packet2);
        
        num_error1to9 = min(n_ber1_1to9);
        num_error1 = min(n_ber1_A2L);
        
%         if num_error2<num_error1to9
%             disp('ERROR! Decoding interference!')
%             disp(['Rx 1 received packet #1to9 error bit number is ',num2str(num_error1to9),'!'])
%             disp(['Rx 1 received ERROR #A2L error bit number is ',num2str(num_error2),'!'])
%             
%         else
%             disp(['Rx 1 received packet #1to9 error bit number is ',num2str(num_error1to9),'!']);
%         end
    else
        disp('Not Detect!');
        num_error1to9 = 55;
    end
    
    [message2] = LoRa_Demodulate_normal(Rx_post_decode2,2);
    if message2 ~= 0
        [n_ber2_A2L] = calculate_ber(message2,packet2);
        [n_ber2_1to9] = calculate_ber(message2,packet1);
        num_errorA2L = min(n_ber2_A2L);
        num_error2 = min(n_ber2_1to9);
%         if num_error2<num_errorA2L
%             disp('ERROR! Decoding interference!')
%             disp(['Rx 2 received packet #A2L error bit number is ',num2str(num_errorA2L),'!']);
%             %         disp(['Rx 2 received ERROR #1to9 error bit number is ',num2str(num_error2),'!']);
%         else
%             disp(['Rx 2 received packet #A2L error bit number is ',num2str(num_errorA2L),'!']);
%         end
    else
        disp('Not Detect!');
        num_errorA2L = 55;
    end
    n_1to9 = min([n_ber1_1to9,num_error2]);
    n_a2l = min([num_errorA2L,num_error1]);
    error_num(pkt_id,:) = [n_1to9,n_a2l];
    disp(['Error num: ',num2str(error_num(pkt_id,1)),' ',num2str(error_num(pkt_id,2))]);
end
filename_str = 'experiment_2x2';
append_to_excel(filename_str, error_num)

end





%-------------------------------------------------------------------------------------%
%%                            SUB- FUNCTION
%-------------------------------------------------------------------------------------%
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
% title('Four-Antenna Beamforming Pattern at 915 MHz with Hamming Window');
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
INPORT_ROOT  = [FileDIR,'input/'];
EXPORT_ROOT  = [FileDIR, 'output/Cluster/'];
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
    int_pkt = int_pkt(1:end-1);
end
for ii = 1 : num_pkt1
    error1 = abs(message(ii+[0:length(int_pkt)-1])-int_pkt);
    n_ber(ii) = length(find(error1>3));
end
ber = n_ber/length(int_pkt);

end

%-------------------------------------------------------------------------------------%
%%
%%
function load_script_flag()
dir = pwd;
Dir = [dir,'/SystemPara/'];
filename = [Dir, 'Script_check_flag.json'];

fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file for reading: %s', filename);
end

% 读取文件内容
raw = fread(fid, inf, 'char=>char');
fclose(fid);

% 将 JSON 字符串解码为结构体
data = jsondecode(raw');

% 提取标志
flag = data.Script_flag;
if flag ==1
    disp('Value == 1!')
elseif flag ==0
    disp('Flag is not set yet...')
end
end
%%
%%
function save_script_flag(flag)
dir = pwd;
Dir = [dir,'/SystemPara/'];
filename = [Dir, 'Script_check_flag.json'];

data.Script_flag = flag;

% 将结构体编码为 JSON 字符串
jsonString = jsonencode(data);

% 将 JSON 字符串写入文件
fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file for writing: %s', filename);
end
fwrite(fid, jsonString, 'char');
fclose(fid);
end

%-------------------------------------------------------------------------------------%
function read_scriptFlag_xml();
%-------------------------------------------------------------------------------------%
%
dir = pwd;
Dir = [dir,'/SystemPara/'];
xmlFileName = [Dir, 'Script_check_flag.xml'];
docNode = xmlread(xmlFileName);

%
flagNode = docNode.getElementsByTagName('Script_flag').item(0);
flag = str2double(flagNode.getFirstChild.getData);
% %
% if flag == 1
%     disp('Flag is set to 1, continuing with processing...');
%     % (假设这里有你的主要代码部分)
% else
%     disp('Flag is not set to 1, aborting...');
% end
xmlFileName = 'Script_check_flag.xml';
xmlwrite([Dir,xmlFileName], docNode);

end
%-------------------------------------------------------------------------------------%

%-------------------------------------------------------------------------------------%
function write_scriptFlag_xml(flag);
%-------------------------------------------------------------------------------------%
%
docNode = com.mathworks.xml.XMLUtils.createDocument('root');
docRootNode = docNode.getDocumentElement;
entryNode = docNode.createElement('Script_flag');
entryNode.appendChild(docNode.createTextNode(num2str(flag)));
docRootNode.appendChild(entryNode);
dir = pwd;
Dir = [dir,'/SystemPara/'];
xmlFileName = 'Script_check_flag.xml';
xmlwrite([Dir,xmlFileName], docNode);
disp('Flag is set to 0, continuing with processing...');
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
        disp('Error: can not find preamble!')
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

disp('Detect preamble, SUCCESS!');
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
        disp('Error: can not find sfd!')
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
        
        disp(['can not find SFD']) ;SFD_offset = 0;
        break;
    end
    st = st + chirplen;
end
SFD_edge = st + SFD_offset;
disp('Detect SFD, SUCCESS!');

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

function extract_signalframe_location(raw_data,EXPORT_ROOT,threshold,DISPLAY_SIGNAL,pkt_len, MAXSearchRange,CoarseStep,FRMlen,rxid);
%==========================================================================================%
% This main function is aim to auto-extract data  (OFDM_LoRa Combined packet)
% OUTPUT:
%       dir: RAWDATA  'USRP_data / Combined_pkt' 
%            /output/
%       frameLocs_cluster1-rx1:  [cluster_information ; relative location]
%                                cluster_information = [cluster_id, cluster_start]
%                                relative location   = cluster_location - cluster_start
%       save as .csv
% check by Ruonan Li  Dec 28 2023

%=================================================================%
%%                   Extract the frame data
%=================================================================%
[SF, BW,Fs] = Parameter_xml_Document_Read();

chirplen = 2 ^ SF;
if nargin == 3
%     threshold      = 3e-2;
    Fs                 = 1e6;
    DISPLAY_SIGNAL     = 1;
    MAXSearchRange     = chirplen* 200;
    pkt_len            = 2 * chirplen ;  % 64 NFFT
    FRMlen             = Fs * 3;
    CoarseStep         = pkt_len*3;
end
if nargin == 4
    Fs =1e6;
    DISPLAY_SIGNAL = 1;
    MAXSearchRange     = chirplen* 50;
    pkt_len            = 45 * chirplen ;  % 64 NFFT
    FRMlen             = Fs * 3;
    CoarseStep         = pkt_len*5;
end
if nargin == 5
   Fs =0.5e6;
   MAXSearchRange     = chirplen* 50;
   pkt_len            = 20 * chirplen ;  % 64 NFFT
   FRMlen             = Fs * 3;
   CoarseStep         = pkt_len*5;
end

% parameter setting

% CoarseStep         = pkt_len*3;
FRM_GAPlen         = FRMlen * 1;
%
PreciseStep        = pkt_len/2;
FRM_LOC_GUARD      = pkt_len/2;
% MAXSearchRange     = chirplen* 50;

% Flag setting
SearchMode         = 1;    % 0: coarse ; 1: precise; -1: nothing to do
FRM_ExtractMode    = 1;    % 0: not extract ; 1: do extract
FRM_CLUSTER_ID     = 1;    % cluster id of frames for signal extracting
% data section

if (SearchMode > 0)
    
    frmLocs_coarse = util_coarse_frameLoc_search(raw_data, CoarseStep,threshold);
    frm_clusters   = util_frame_cluster_extract(frmLocs_coarse, FRM_GAPlen);
    %% select a cluster
    if (SearchMode >= 1)
        cluster_id  = FRM_CLUSTER_ID;
        
        frmLocs_precise = util_precise_frameLoc_extract(...
            cluster_id, frm_clusters, frmLocs_coarse, ...
            raw_data, MAXSearchRange, FRM_LOC_GUARD, threshold, PreciseStep);
        % extract the signal of frames in target cluster,
        % write into separated files
        if (FRM_ExtractMode > 0)
            util_frameCluster_signal_extract(raw_data, frmLocs_precise, cluster_id,EXPORT_ROOT,rxid);
        end
    else
        frmLocs_precise = [frmLocs_coarse];
    end
    if (DISPLAY_SIGNAL > 0)
        %=================================================================%
        figure;
        plot(abs(raw_data),'k'); hold on
        % draw frame locations
        sz      = size(frmLocs_precise);
        num_frm = sz(1);
        
        y       = max(abs(raw_data));
        for idx = 1:num_frm
            
            x1  = frmLocs_precise(idx, 1);
            plot([x1 x1], [0 y], '-r', 'linewidth', 1);
            
            x2  = frmLocs_precise(idx, 2);
            plot([x2 x2], [0 y], '-r', 'linewidth', 1);
        end
        
        xlabel('Samples');
        ylabel('Amplitude');
        title('Rx Raw signals')
        %=================================================================%        
    end
end
end

%----------------------------------------------------------------------------------------------------------%
%%                                 utility functions
%----------------------------------------------------------------------------------------------------------%
%----------------------------------------------------------------------------------------------------------%
function util_frameCluster_signal_extract(sigs, cluster_frameLocs, cluster_id,DIR,rxid)
%   this function extract all frame signal and get their start and end location.
%   INPUT:
%        sigs:               raw signal
%        cluster_frameLocs:  each cluster start and end location [start end]
%        cluster_id:
%        RX_ID
%        DIR:                file save root.
%   OUTPUT:
%        frameSigs_cluster   only each extracted frame signal
%        frameLocs_cluster   [cluster_info; relative_locs];
%                            cluster_info  = [cluster_id, first cluster start bin];
%                            relative_locs = [cluster_frameLocs - cluster_st].
%
sz      = size(cluster_frameLocs);
num_frm = sz(1);
cluster_st = cluster_frameLocs(1, 1);
cluster_ed = cluster_frameLocs(num_frm, 2);
if (cluster_ed > numel(sigs)) 
    
    cluster_ed = numel(sigs);
end
sig_len           = cluster_ed - cluster_st + 1;

cluster_sig       = sigs(cluster_st:cluster_ed);
sig_data          = zeros(2*sig_len, 1);
sig_data(1:2:end) = real(cluster_sig);
sig_data(2:2:end) = imag(cluster_sig);

% record signal info
cluster_info  = [cluster_id, cluster_st-1];
relative_locs = cluster_frameLocs - cluster_st+1;
sig_info      = [cluster_info; relative_locs];

file_info     = [DIR 'frameLocs_cluster_Rx' num2str(rxid) , '.csv'];
dlmwrite(file_info, sig_info, 'precision', '%+.20e');
end
%----------------------------------------------------------------------------------------------------------%

%----------------------------------------------------------------------------------------------------------%
function frame_locs = util_precise_frameLoc_extract(...
    cluster_id, clusters, coarse_frameLocs, ...
    sigs, frm_max_len, frm_guard_len, pwr_threshold, mov_win)
%   this function is to get the frame signal start and end location after
%   precise search
%   INPUT:
%         cluster_id:
%         clusters:         cluster number [1 x]
%         coarse_frameLocs: after coarse search the location might be..
%         sigs:             raw signal
%         frm_max_len:      maximum search length
%         frm_guard_len:    guard window length
%         pwr_threshold:    true value
%         mov_win:          precise step
%   OUTPUT:
%         frame_locs:  each cluster start and end location [start end]

num_sig     = numel(sigs);
sz          = size(clusters);
num_cluster = sz(1);
if (cluster_id <= 0 || cluster_id > num_cluster)
    frame_locs = [];
    return;
end

cluster_start_idx = clusters(cluster_id, 1);
cluster_frame_num = clusters(cluster_id, 2);

frame_locs = zeros(cluster_frame_num, 2);   % start and end pos of frames in the raw signals
for idx = 1:cluster_frame_num
    frame_idx  = cluster_start_idx + idx - 1;
    
    % search for precise frame locations
    coarse_loc = coarse_frameLocs(frame_idx);
    [frm_loc, frm_len] = util_beacon_locate(...
        sigs, coarse_loc, frm_max_len, pwr_threshold, mov_win);
    
    frm_loc_st = frm_loc - frm_guard_len;
    frm_loc_ed = frm_loc + frm_len + frm_guard_len;
    if (frm_loc_st <= 0)
        frm_loc_st = 1;
    end
    if (frm_loc_ed > num_sig)
        frm_loc_ed = num_sig;
    end
    
    frame_locs(idx, :) = [frm_loc_st, frm_loc_ed];
end
end
%----------------------------------------------------------------------------------------------------------%

%----------------------------------------------------------------------------------------------------------%
function clusters = util_frame_cluster_extract(coarse_frameLocs, cluster_gap)
% this function is to extract the cluster after coarse search
%       INPUT:
%           coarse_frameLocs:
%           cluster_gap:
%       OUTPUT:
% clusters = [start frame index, # of frames in cluster]

num_frames = numel(coarse_frameLocs);

if (num_frames <= 0)
    clusters = [];
    return;
end


clusters    = zeros(num_frames, 2);
num_cluster = 0;

cluster_start_index = 1;    % the start frame of the current cluster
cluster_frame_num = 1;      % # of frames in the current cluster
for idx = 2:num_frames
    if (coarse_frameLocs(idx) - coarse_frameLocs(idx-1) > cluster_gap)
        % a new cluster detected
        num_cluster = num_cluster + 1;
        clusters(num_cluster, :) = [cluster_start_index, cluster_frame_num];
        
        % initialize parameters for the new cluster
        cluster_start_index = idx;
        cluster_frame_num   = 1;
    else
        cluster_frame_num   = cluster_frame_num + 1;
    end
end

% handle the last cluster
num_cluster = num_cluster + 1;
clusters(num_cluster, :) = [cluster_start_index, cluster_frame_num];

clusters    = clusters(1:num_cluster, :);
end
%----------------------------------------------------------------------------------------------------------%

%----------------------------------------------------------------------------------------------------------%
function sig_locs = util_coarse_frameLoc_search(sigs, step_win, pwr_threshold)
% this function is to coarse search
%       INPUT:
%           sigs:           raw data
%           step_win:       coarse step length
%           pwr_threshold:  true value
%       OUTPUT:
%           sig_locs:       after coarse search the frame location might be..

sigs        = abs(sigs);
num_samples = numel(sigs);
num_win     = ceil(num_samples/step_win);

% indicate the presence of chirp signal in a specific win
sigIndicator_wins = zeros(num_win, 1);

for pos = 1:step_win:num_samples
    win_st = pos;
    win_ed = pos + step_win-1;
    if (win_ed > num_samples)
        win_ed = num_samples;
    end
    
    if (mean(sigs(win_st:win_ed)) > pwr_threshold)
        win_id = floor(win_st/step_win) + 1;
        sigIndicator_wins(win_id) = 1;
    end
end

% detect the coarse location of frames
frm_stLocs = zeros(num_win, 1);
frm_num    = 0;

for w = 1:num_win-1
    if ((w == 1 && sigIndicator_wins(w) > 0) || ...
            (sigIndicator_wins(w) == 0 && sigIndicator_wins(w+1) > 0))
        % presence of signal detected
        frm_num = frm_num + 1;
        
        % estimate the starting signal position
        sig_pos = (w-2) * step_win;
        if (sig_pos < 1)
            sig_pos = 1;
        end
        
        frm_stLocs(frm_num) = sig_pos;
    end
end

sig_locs = frm_stLocs(1:frm_num, :);
end
%----------------------------------------------------------------------------------------------------------%

%----------------------------------------------------------------------------------------------------------%
function [frm_st, frm_len] = util_beacon_locate(...
    signals, init_pos, search_len, pwr_threshold, mov_win)
%  find the precise start bin
%       INPUT:
%             signals:       raw signal
%             init_pos:      after coarse search, the might frame location
%             search_len:    coarse search range
%             pwr_threshold: true value
%             mov_win:       precise step
%
%       OUTPUT:
%             frm_st:        frame start bin
%             frm_len:       each frame length

sig_num = numel(signals);
ed_pos  = init_pos + search_len;
if (ed_pos > sig_num)
    ed_pos = sig_num;
end

win_half = ceil(mov_win/2);
mag      = abs(signals(init_pos:ed_pos));
mag      = movmean(mag, win_half);
mag_num  = numel(mag);

frm_st  = 0;
frm_len = 0;
for pos = 1:mag_num
    if (mag(pos) > pwr_threshold)
        if (frm_st  <= 0)
            frm_st  = pos;
            frm_len = 1;
        else
            frm_len = pos - frm_st + 1;
        end
    elseif (frm_len > 0) % the first frame ends
        break;
    end
end

frm_st = frm_st + init_pos - 1;
end
%----------------------------------------------------------------------------------------------------------%

%----------------------------------------------------------------------------------------------------------%
function res = movmean(array, range)
% this function is to do movement mean() in a range
%       INPUT:
%          array: extracted data
%          range:
%       OUTPUT:
%          res: move window mean result

len = numel(array);

res = zeros(size(array));
for idx = 1:len
    st  = idx - range;
    ed  = idx + range;
    if (st <= 0)
        st = 1;
    end
    if (ed > len)
        ed = len;
    end
    
    res(idx) = mean(array(st:ed));
end
end
%----------------------------------------------------------------------------------------------------------%



