function exp_eval_plot_Spatial_Reuse_2_tx_CDF_PRR()
% for the Sec2.2 deploy MIMO_LoRa
clc;% close all
SF = 8;
BW = 125e3;
Fs = BW;
%----------------------------------W / O CFO-------------------------------%

% ===================== %
%   Read cots node .txt
% ===================== %

num_packet = 50;
folder_vary_a2l = [1:4];
folder_vary_1to9 = [6:10];
folder_vary=[1:4,6:10];
%----------------------------------W -------------------------------%

for ii = 1 : 10
    num = ii;
    for folder_id = 1:9
        str = ['m', num2str(folder_vary(folder_id))] ;
        foldername = folder_vary(folder_id);

        [COTS_result,num_pkt(ii,folder_id)]    = Read_COTS_saved_TXT(str,num);
        if ~isempty(COTS_result)
        [rate_a2l,rate_1to9] = COTS_correct_string(COTS_result);
        correct_idx_1to9     = find(rate_1to9< 0.78);
        correct_num(folder_id,1) = length(correct_idx_1to9);
        correct_idx_a2l      = find(rate_a2l< 0.78);
        correct_num(folder_id,2)  = length(correct_idx_a2l);
        if foldername == 1||foldername == 2||foldername == 3||foldername == 9||foldername == 10
            prr_num(folder_id) =  correct_num(folder_id,1) ; % 1to9
        elseif foldername == 4||foldername == 6||foldername == 7||foldername == 8
            prr_num(folder_id) =  correct_num(folder_id,2) ; % a2l
        end
        else
            prr_num(folder_id) = 0;
        end
    end
    prr(:,ii) = prr_num/min(num_pkt(ii,:));
end
%-----------------------------  Read 11 12 (20min)  ---------------------------------%
n_group = 10;
group_len = 30;
for ii = 11 : 12
    num = ii;
    for folder_id = 1:9
        str = ['m', num2str(folder_vary(folder_id))] ;
        foldername = folder_vary(folder_id);

        [COTS_result,num_pkt(ii,folder_id)]    = Read_COTS_saved_TXT(str,num);
        if ~isempty(COTS_result)
            COTS_result_select = {};
            for groupid = 1:n_group
                if length(COTS_result)<10*group_len*2
%                     COTS_result = {COTS_result;COTS_result{20*group_len-length(COTS_result)}};
                    group_len = 20;
                end
                for pkt = 1:group_len*2
                    COTS_result_select{pkt,1} = [COTS_result{(groupid-1)*group_len*2+pkt,1}];
                end

                [rate_a2l,rate_1to9] = COTS_correct_string(COTS_result_select);
                correct_idx_1to9     = find(rate_1to9< 0.78);
                correct_long(groupid,1) = length(correct_idx_1to9);
                correct_idx_a2l      = find(rate_a2l< 0.78);
                correct_long(groupid,2)  = length(correct_idx_a2l);
                if foldername == 2||foldername == 3||foldername == 9||foldername == 10
                    prr_long(groupid,:) =  correct_long(groupid,1)/group_len; % 1to9
                elseif foldername == 1||foldername == 4||foldername == 6||foldername == 7||foldername == 8
                    prr_long(groupid,:) =  correct_long(groupid,2)/group_len ; % a2l
                end
            end
        else
            prr_long(folder_id) = 0;
        end
        prr(folder_id,(ii-10)*10+[1:10]) = prr_long;

    end
        
end
tx_2_prr = prr;
tx_2_prr(find(prr>1)) = 1;





save('2tx_prr_9_30times','tx_2_prr');

end
%-------------------------------------------------------------------------------------%
%%                        Sub - Function
%-------------------------------------------------------------------------------------%
% function [data] = OpenSource_data_biterrornum(time_para)
% FileHead  = [time_para,'_BER_opensource_50times','.xlsx'];
% dir       = pwd;
% Dir_all   = [dir,'/',time_para,'_COTSnode/'];
% file_info = [Dir_all,FileHead];
% data     = readtable(file_info);
% disp(['Read ', FileHead, ' data...' ])
% end
%-------------------------------------------------------------------------------------%
%-------------------------------------------------------------------------------------%
function [charDATA,num_pkt] = Read_COTS_saved_TXT(str,num)
%str: 1_1to9; 2_A2L
SubFolderHead= [str,'\2tx'];
dir       = pwd;
Dir_all   = [dir,'\Spatial_reuse_opportunity_RawData\Spatial_reuse_opportunity\',SubFolderHead,'\'];
file_info = [Dir_all,num2str(num),'.txt'];
fi_1     = fopen(file_info,'r','n', 'UTF-8');
charDATA = cell({});
num_pkt = 0;
while 1
    values   = fgetl(fi_1);
    if ~ischar(values)
        break;
    end
%     if values()
    charDATA  =[charDATA; values];
  
end
  % 检查文件是否存在
   
fclose(fi_1);

if isempty(charDATA)
    disp([SubFolderHead, 'Not detected the Message!'])
else
    disp(['Read ', SubFolderHead,num2str(num), ' data...' ])
end
for ii = 1:length(charDATA)
    value1 = charDATA{ii};
    if contains(value1 ,'SER')
        num_pkt = num_pkt+1;
    end
end
end
%-------------------------------------------------------------------------------------%

function [rate_a2l,rate_1to9] = COTS_correct_string(data_text)

Str1  = 'ABCDEFGHJKL';
Str2  = '123456789#';
Header= 'SERVER got packet from client: ';
comparedStr1   = [Header,Str1];
comparedStr2   = [Header,Str2];
lockon         = strfind(data_text, Header);
lockon_num     = length(find(cellfun('length',lockon) ~=0));
for rr = 1 : size(data_text,1)
    data_time = char(data_text{rr});
    if data_time(1)~='S'
    [rate_a2l(rr,:)]  = mystrsim(comparedStr1,data_time(17:end));
    [rate_1to9(rr,:)] = mystrsim(comparedStr2,data_time(17:end));
    else
    [rate_a2l(rr,:)]  = mystrsim(comparedStr1,data_time(1:end));
    [rate_1to9(rr,:)] = mystrsim(comparedStr2,data_time(1:end));
    end
end
end
%-------------------------------------------------------------------------------------%

function strsim = mystrsim(target,source)
len1= length(source);
len2=length(target);
d=zeros(len1+1,len2+1);
for i=2:length(source)+1
    d(i,1)= i;
end
for j=2:length(target)+1
    d(1,j)= j;
end
for i =2:length(source)+1
    for j = 2:length(target)+1
        if source(i - 1) == target(j - 1)
            d(i,j) = d(i-1,j-1);
        else
            edIns = d(i,j-1)+1;
            edDel = d(i-1,j)+1;
            edRep = d(i-1,j-1)+1;
            d(i,j)= min(min(edIns, edDel),edRep);
        end
    end
end
y=d(length(source)+1,length(target)+1);
strsim=1/(exp(1./y));

end

%-------------------------------------------------------------------------------------%
function  fig_data_prr_errorbar(prr_data)

% 计算每组数据的均值和标准误差
mean_prr = mean(prr_data, 3);
std_prr = std(prr_data, 0, 3);
num_packets = size(prr_data, 3);
stderr_prr = std_prr / sqrt(num_packets);

% 设置组名称和图例名称
groupsName = {'Group 1', 'Group 2', 'Group 3'};
legendName = {'Link 1', 'Link 2'};

% 创建分组柱状图
figure;
b = bar(mean_prr, 'grouped');
hold on;

% 设置柱的颜色为灰度
b(1).FaceColor = [0.2, 0.2, 0.2]; % 灰度值
b(2).FaceColor = [0.6, 0.6, 0.6]; % 灰度值

% 添加误差条
ngroups = size(groupsName, 2);
nbars = size(legendName, 2);

% 计算每个组的位置
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    % 计算每个柱的位置
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, mean_prr(:, i), stderr_prr(:, i), 'k', 'linestyle', 'none');
end

% 添加图例
legend(legendName, 'FontName', 'Arial', 'fontsize', 12, 'location', 'northeast');

% 设置x轴标签和其他属性
set(gca, 'xticklabel', groupsName, 'box', 'off', 'FontName', 'Arial', 'fontsize', 11, 'xtick', 1:ngroups);
title('Packet Reception Ratio with Error Bars');
xlabel('Groups', 'fontsize', 12);
ylabel('PRR', 'fontsize', 12);
hold off;









end