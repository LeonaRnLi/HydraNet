function exp_sec5_RSS_calculate()
 close all

%% ========  Pre SET Parameters  ======= %
% Read EVM file
dir          = pwd;
SF           = 8;
Fs           = 125e3;
n_symbol     = 2^SF;
%
%% ============= STEP 1 ================ %
%  Extract the packet
%  we will not estimate the first packet.
%  we can choose the second or third.
BW = Fs;Vref = 1/2^14;
[raw_data,EXPORT_ROOT] = GetDir_Read_lora_Data;
% figure;plot(abs(raw_data));
threshold      = 3e-2;
extract_signalframe_location(raw_data,EXPORT_ROOT, threshold);
cluster_loc    = ['frameLocs_cluster1', '.csv'];
cluster_info   = csvread([EXPORT_ROOT cluster_loc]);
startnum       = cluster_info(1,2);
num_frame      = size(cluster_info,1)-1;
for pkt = 1 : num_frame
    relativeID     = pkt+1;
    cluster_st     = startnum + cluster_info(relativeID,1)+100;
    cluster_ed     = cluster_st+n_symbol * 3;
    rx_signal_detect = raw_data(cluster_st:cluster_ed);
%         FigTFSpectrum(rx_signal_detect,Fs);
    n_chirp   = floor(length(rx_signal_detect)/n_symbol);
    up_chirp  = loramod(zeros(1,n_chirp),SF,BW,Fs,1) ; % upchirp
    buffer    = up_chirp.*rx_signal_detect(1:n_chirp*n_symbol);
    down_edge = fftshift(fft(reshape(buffer,n_symbol,n_chirp)),1);
%     figure;imagesc(abs(down_edge))
    max_val   = max(max(down_edge));
    norm_down = down_edge/max_val;
    [a index0] = find(abs(norm_down)>0.96);
    index = unique(index0);
    
    sig    = rx_signal_detect(1:2*n_symbol+n_symbol);
    %    find start edge -- correlation
    
    %         FigTFSpectrum(sig,Fs);
    [edges,val] = trace_preamble_edge_correlation(sig,SF,BW,Fs,1,-1);
    temp_pream = edges(find(val>=0.4));
    LoRa_edge_idx(pkt) = index0(1)+temp_pream(1)-1;
    temp = LoRa_edge_idx(pkt);n_len = 2;
        %         FigTFSpectrum(rx_signal_detect(LoRa_edge_idx(pkt,1)+[0:256]),Fs);
     for cc = 1:1
         sig = rx_signal_detect(temp+[0:256*n_len]);
     LoRa_RSS(pkt,cc) = mean(10*log10(abs(sig*Vref).^2));
%      FigTFSpectrum(rx_signal_detect(temp+[0:256]),Fs);
     temp = temp+256*n_len;
     end
  disp(['RSS: ', num2str(LoRa_RSS(pkt,:)),'dBm'])
  disp(['Link fading: ', num2str(23-LoRa_RSS(pkt,:)),'dB!'])
end





end
%-------------------------------------------------------------------------------------%
%%      Sub - function  
%-------------------------------------------------------------------------------------%
function [raw_data,EXPORT_ROOT] = GetDir_Read_lora_Data();
FileHead = 'USRP_RSS_TEST';
Filestr = [FileHead];
DataDIR  = '/RAWDATA/USRP_RSS/';
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
EXPORT_ROOT = [dir,'/RAWDATA/USRP_RSS/','output/'];
end
%-------------------------------------------------------------------------------------%
%-------------------------------------------------------------------------------------%

function [edges,val] = trace_preamble_edge_correlation(sig,SF,BW,Fs,flag,display_flag);
%   flag > 0  find preamble edge
%   flag < 0  find sfd edge
n_symbol = 2^SF* Fs/BW;
% [up_chirp, down_chirp] = chirp_template(n_symbol);
if flag > 0
        base_chirp =  loramod(0,SF,BW,Fs,1) ; % downchirp
%     base_chirp = down_chirp(1:n_symbol);
elseif flag < 0
        base_chirp =  loramod(0,SF,BW,Fs,-1) ; % upchirp
%     base_chirp = up_chirp(1:n_symbol);
end
    fft_loc_num = numel(sig) - n_symbol + 1;
    corr_res = -1*ones(fft_loc_num, 1);
     for loc = 1:fft_loc_num
        seg_st = loc;
        seg_ed = seg_st + n_symbol - 1;
        seg = sig(seg_st:seg_ed);
%         FigTFSpectrum(seg,Fs); 
        R = corrcoef(seg, base_chirp);
        corr_res(loc) = abs(R(1,2));
%      
    end
      % detect correlation peaks
    edges = zeros(ceil(fft_loc_num/n_symbol), 1);
    val   = zeros(size(edges));
    edge_num = 0;
    for idx = 1:fft_loc_num
        if (corr_res(idx) > 0.2)
            edge_num = edge_num + 1;
            edges(edge_num) = idx;
            val(edge_num) = corr_res(idx);
        end
    end
    edges = edges(1:edge_num);
    val   = val(1:edge_num);
    if display_flag > 0
    fig_preamble_correlation(corr_res);
    end


end
%-------------------------------------------------------------------------------------%
%-------------------------------------------------------------------------------------%

function fig_preamble_correlation(corr_res)

    Y = corr_res;
    
    figure
    plot(Y, 'color', 0.2*[1 1 1],  'linewidth', 1);
%     stem(Y, '.', 'color', 0.2*[1 1 1], 'linewidth', 2);
    
    fontsize = 14;
    fig_width = 16;
    fig_height = 5;
    xlabel('PHY sample #');
    %xlabel('Time (\mus)');
    ylabel('Normalized correlation');
    title('');
    
    set(gca, 'fontsize', fontsize, 'fontname', 'Times New Roman');

    set(gcf,'unit', 'centimeters', 'position', [5 5 fig_width fig_height]);
    left_margin = 0.10;
    right_margin = 0.05;
    bot_margin = 0.23; 
    top_margin = 0.1;
    set(gca, 'position', [left_margin, bot_margin, 1-left_margin-right_margin, 1-bot_margin-top_margin]);
    

    X_max = numel(corr_res);
    xlim([0 X_max]);
    ylim([0 1]);
    XTicks = floor([1 (0.2:0.2:1)*X_max]);
    XLabels = XTicks; %* 1/250 * 1e3;
    set(gca, 'XTick', XTicks, 'XTickLabel', XLabels, 'fontsize', fontsize-2);
    set(gca, 'YTick', 0:0.2:1, 'YTickLabel', {'0.0' '0.2' '0.4' '0.6' '0.8' '1.0'}, 'fontsize', fontsize-2);
    
    grid on
    set(gca,'gridlinestyle', ':', 'gridcolor', 0.0*[1 1 1]);
    
    datacursormode on
end
%-------------------------------------------------------------------------------------%



%-------------------------------------------------------------------------------------%

function [y] = loramod(x,SF,BW,fs,varargin)
% loramod LoRa modulates a symbol vector specified by x
%
%   in:  x          1xN symbol vector wher N=1-Inf
%                   with values {0,1,2,...,2^(SF)-1}
%        BW         signal bandwidth of LoRa transmisson
%        SF         spreading factor
%        Fs         sampling frequency
%        varargin{1} polarity of chirp
%
%  out:  y          LoRa IQ waveform
if (nargin < 4)
    error(message('comm:pskmod:numarg1'));
end

if (nargin > 5)
    error(message('comm:pskmod:numarg2'));
end

% Check that x is a positive integer
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

gamma   = x/Ts ;     % f0
beta    = BW/Ts ; %  Kr

time    = (0:Ns-1)'.*1/fs ;
freq    = mod(gamma + Inv.*beta.*time,BW) - BW/2 ;
Theta   = cumtrapz(time,freq) ;
% figure; plot(freq(:,2)/pi)
y       = reshape(exp(j.*2.*pi.*Theta),numel(Theta),1) ;
end
%-------------------------------------------------------------------------------------%
%-------------------------------------------------------------------------------------------------------%
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

figure(1);
overlap   = floor(win*7/8);
[S,F,T,P] = spectrogram(signal, win, overlap, nfft, Fs, 'centered');
norm_freq = F/1000;     % unit: kHz
win_id    = 1:1:numel(T);  % window ID
if (LOG_PSD == 1)
    psd = 10*log10(abs(P));
else
    psd = abs(P);
end
fig = surf(win_id, norm_freq, psd);
axis tight; shading interp; view(0, 90);
ylim([-1 1] * max(norm_freq));
xlabel('Win ID'); ylabel('Frequency (kHz)');
title(str); grid off;
%colorbar
set(gcf, 'unit', 'centimeters', 'position', [1/4, 8, 40, 5]);
left_margin  = 0.05;
right_margin = 0.02;
bot_margin   = 0.22;
top_margin   = 0.15;
set(gca, 'position', [left_margin, bot_margin, 1-left_margin-right_margin, 1-bot_margin-top_margin]);

end
%-------------------------------------------------------------------------------------%
