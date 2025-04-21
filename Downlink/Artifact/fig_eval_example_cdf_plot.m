% Example data
% Define file paths
% file_list = {'Baseline_prr_9_30times.mat', '2tx_prr_9_30times.mat', '3tx_prr_9_30times.mat', '4tx_prr_9_30times.mat'}; % Replace with actual file paths
% close all
% Initialize PRR_matrix
PRR_matrix = zeros(4, 30);

% Compute average PRR
for i = 1:length(file_list)
    data = load(file_list{i});
    % Get variable name in the .mat file
    var_info = whos('-file', file_list{i});
    var_name = var_info.name;
    matrix = data.(var_name);
    % Compute average of each column
    PRR_matrix(i, :) = mean(matrix, 1);
end

payload_size = 10; % Set payload size to 10 bytes
transmission_time = 0.08; % Set transmission time to 80 milliseconds
filename = 'CDF_throughput_exp5.mat'; % Output file name

% Call function to generate MAT file
generate_throughput_mat(PRR_matrix, payload_size, transmission_time, filename);

data = load('CDF_throughput_exp5.mat');

% Extract data
link1 = data.link1;
link2 = data.link2;
link3 = data.link3;
link4 = data.link4;
links = [data.link1, data.link1]; % LMAC & TPC

% Set legend names
legendNames = {['HydraNet',newline,' (1 Ant.)'], ['HydraNet',newline,' (2 Ant.)'],...
    ['HydraNet',newline,' (3 Ant.)'], ['HydraNet',newline,' (4 Ant.)'],...
  };

% Draw CDF plot
figure;
hold on;

% Predefined colors and line styles
colors = [0.7 0.7 0.7; 0.5 0.5 0.5; 0.3 0.3 0.3; 0.1 0.1 0.1;...
    [177 34 34]/255; [208, 122, 122]/255]; % Different colors for different configs

lineStyles = {'-', '--', '-.', ':', '-', ':'}; % Different line styles
markers = {'o', 's', 'd', 'p', '^', 'v'}; % Different markers
markerFrequency = 4;
handles = [];

for i = 1:4
    link_data = data.(sprintf('link%d', i));
    [f, x] = ecdf(link_data); % Calculate CDF
    markerIndices = 1:markerFrequency:length(x);
    % Plot CDF curve
    p = plot(x, f, 'LineWidth', 2.5, 'Color', colors(i,:), 'LineStyle', lineStyles{i}, ...
        'Marker', markers{i}, 'MarkerIndices', markerIndices, 'DisplayName', legendNames{i}, ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:), 'MarkerSize', 5);
    handles(end+1) = p;
end

% Set figure properties
fontsize = 14;
figureWidth = 9;
figureHeight = 6.7;
customLegendOrder = [1,2,3,4]; % Adjust legend order as needed
customLegendNames = legendNames(customLegendOrder);
lgd = legend(handles(customLegendOrder(1:end)), customLegendNames(1:end),...
    'Location', 'bestoutside', 'NumColumnsMode','manual','NumColumns',1, ...
    'FontName', 'Times New Roman', 'fontsize', fontsize-2, 'Box', 'off'); hold on
lgd.ItemTokenSize = [18,10];

set(gca, 'fontsize', fontsize, 'fontname', 'Times New Roman');
set(gcf, 'unit', 'centimeters', 'position', [20 5 figureWidth figureHeight]);
leftMargin = 0.16;
rightMargin = 0.32;
botMargin = 0.21;
topMargin = 0.07;
set(lgd,...
    'Position', [0.684368468859177 0.0936254690422548 0.320588235294118 0.903162029185314]);

set(gca, 'position', [leftMargin, botMargin, 1-leftMargin-rightMargin, 1-botMargin-topMargin]);
xlabel('Throughput (kbps)', 'FontName', 'Times New Roman', 'fontsize', fontsize);
ylabel('CDF', 'FontName', 'Times New Roman', 'fontsize', fontsize);
ax = gca;
ax.YGrid = 'on';
ax.XGrid = 'on';
ax.GridColor = [104, 109, 118] / 256;
ylim([0 1]);
xlim([0 5]);
XTicks = [0:1:5];
XLabels = arrayfun(@(x) num2str(x), XTicks, 'UniformOutput', false);
set(gca, 'XTick', XTicks, 'XTickLabel', XLabels, 'fontsize', fontsize, 'XTickLabelRotation', 0);
YTicks = 0:0.2:1;
set(gca, 'YTick', YTicks, 'fontsize', fontsize, 'XTickLabelRotation', 0);

ax.LineWidth = 1;
set(gca, 'LineWidth', 2);

box on;
hold off;

dir = pwd;
set(gcf, 'renderer', 'Painters');
print('-depsc2', fullfile(dir, 'Sec5_SpatialReuse_throughput'), '-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function generate_throughput_mat(PRR_matrix, payload_size, transmission_time, filename)
    % Input parameters:
    % PRR_matrix - A 4-row matrix of PRR values (column count can vary)
    % payload_size - Payload size per packet in bytes
    % transmission_time - Transmission time per packet in seconds
    % filename - Output filename for the MAT file

    % Ensure PRR_matrix has 4 rows
    assert(size(PRR_matrix, 1) == 4, 'PRR matrix must have 4 rows');

    % Calculate throughput with configuration-specific scaling factors
    throughput_matrix = zeros(size(PRR_matrix));
    throughput_matrix(1, :) = PRR_matrix(1, :) * 2 * payload_size * 8 / transmission_time / 1000; % 1Tx 
    throughput_matrix(2, :) = PRR_matrix(2, :) * 3.5 * payload_size * 8 / transmission_time / 1000; % 2Tx 
    throughput_matrix(3, :) = PRR_matrix(3, :) * 4.5 * payload_size * 8 / transmission_time / 1000; % 3Tx 
    throughput_matrix(4, :) = PRR_matrix(4, :) * 5.3 * payload_size * 8 / transmission_time / 1000; % 4Tx 

    % Save data to MAT file
    mat_data = struct();
    for i = 1:4
        mat_data.(sprintf('link%d', i)) = throughput_matrix(i, :)';
    end
    save(filename, '-struct', 'mat_data');
    
    % Print confirmation message
    fprintf('Data successfully saved to %s\n', filename);
end
