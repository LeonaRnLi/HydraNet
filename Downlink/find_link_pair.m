function [optimized_Tx_power,optimized_di_theta] = find_link_pair(link_attenuation,link_angle_deg)
clc;
close all
% find link Mi>=Nj gateway>=user to get the optimized weight
% link_attenuation:a11= gateway1 to user1



theta_vary   = linspace(-90, 90, 360); 
num_gateways = size(link_attenuation, 1);
num_users    = size(link_attenuation, 2);
Na           = 2;  % can change  

% --------------------------------- %
%           optimization 
d_theta = -60*ones(1,num_users);
lb      = [1+zeros(1, num_users),d_theta];
ub      = [20 * ones(1, num_users),-d_theta];
options = optimoptions('particleswarm', 'Display', 'none', 'SwarmSize', 100, 'MaxIterations', 400);

[optimized_params, fval] = particleswarm(@(params) -objective_function(params, link_attenuation, link_angle_deg, theta_vary, Na), ...
    2 * num_gateways, lb, ub, options);

optimized_Tx_power = optimized_params(1:num_gateways);
optimized_di_theta = optimized_params(num_gateways+1:end);

% disp('Optimized Tx_power:'); disp(optimized_Tx_power);
% disp('Optimized divation theta:'); disp(optimized_di_theta);

% RSS = calculate_rss(optimized_Tx_power, optimized_di_theta, link_attenuation, link_angle_deg, theta_vary, Na);
% disp('Optimized RSS:'); disp(RSS);

% 计算最终的Cij矩阵
Cij = calculate_cij(optimized_Tx_power, optimized_di_theta, link_attenuation, link_angle_deg, theta_vary, Na);
disp('Selected Gateway and user index matrix Cij:'); disp(Cij);

% --------------------------------- %
%       Output BF weight
angle_deg_real = sum(link_angle_deg.*Cij,2);
theta_real = [optimized_di_theta.'+angle_deg_real];

end

% ------------------------------------------------------------------ %

function score = objective_function(params, link_attenuation, link_angle_deg, theta, Na)
[num_gateways,num_users] = size(link_attenuation);

Tx_power = params(1 : num_gateways);
div_theta = params(num_gateways+1 : end); 
RSS = calculate_rss(Tx_power,div_theta, link_attenuation, link_angle_deg, theta, Na);

Cij = zeros(num_gateways, num_users);

% 计算Cij矩阵和约束条件
for user = 1:num_users
    best_gateway = -1;
    best_rss = -Inf;
    for gw = 1:num_gateways
        if Tx_power(gw) >= 0 && Tx_power(gw) <= 20 % 发射功率约束
            if RSS(gw,user)>=-120
                other_gateways = setdiff(1:num_gateways, gw);
                interference = (RSS(other_gateways, user));
                if RSS(gw, user) > (interference + 0.8) % 干扰约束
                    if RSS(gw, user) > best_rss
                        best_rss = RSS(gw, user);
                        best_gateway = gw;
                    end
                end
            end
        end
    end
    if best_gateway > 0
        Cij(best_gateway, user) = 1;
    end
end

% 保证Cij矩阵每个网关只能与一个用户通信
for gw = 1:num_gateways
    if sum(Cij(gw, :)) > 1
        [~, idx] = max(RSS(gw, :));
        Cij(gw, :) = 0;
        Cij(gw, idx) = 1;
    end
end

% 保证Cij矩阵每行每列均只有一个1
for user = 1:num_users
    if sum(Cij(:, user)) > 1
        [~, idx] = max(RSS(:, user));
        Cij(:, user) = 0;
        Cij(idx, user) = 1;
    end
end

% 目标函数：最大化Cij矩阵中1的数量
score = sum(Cij(:));

end
% ------------------------------------------------------------------ %
function Cij = calculate_cij(Tx_power, div_theta, link_attenuation, link_angle_deg, theta, Na)
[num_gateways, num_users] = size(link_attenuation);

RSS = calculate_rss(Tx_power, div_theta, link_attenuation, link_angle_deg, theta, Na);
Cij = zeros(num_gateways, num_users);

% 计算Cij矩阵和约束条件
for user = 1:num_users
    best_gateway = -1;
    best_rss = -Inf;
    for gw = 1:num_gateways
        if Tx_power(gw) >= 0 && Tx_power(gw) <= 20 % 发射功率约束
            if RSS(gw,user)>=-120
                other_gateways = setdiff(1:num_gateways, gw);
                interference = RSS(other_gateways, user);
                if RSS(gw, user) > (interference + 0.8) % 干扰约束
                    if RSS(gw, user) > best_rss
                        best_rss = RSS(gw, user);
                        best_gateway = gw;
                    end
                end
            end
        end
    end
    if best_gateway > 0
        Cij(best_gateway, user) = 1;
    end
end

% 保证Cij矩阵每个网关只能与一个用户通信
for gw = 1:num_gateways
    if sum(Cij(gw, :)) > 1
        [~, idx] = max(RSS(gw, :));
        Cij(gw, :) = 0;
        Cij(gw, idx) = 1;
    end
end

% 保证Cij矩阵每行每列均只有一个1
for user = 1:num_users
    if sum(Cij(:, user)) > 1
        [~, idx] = max(RSS(:, user));
        Cij(:, user) = 0;
        Cij(idx, user) = 1;
    end
end

end
% ------------------------------------------------------------------ %

% ------------------------------------------------------------------ %

function RSS = calculate_rss(Tx_power,div_theta, link_attenuation, link_angle_deg, theta, N)
[num_gateways,num_users] = size(link_attenuation);
div_theta = div_theta.';
% initial the RSS matrix
RSS = zeros(num_gateways, num_users);

for gw = 1:num_gateways
    for user = 1:num_users
        [~, idx] = min(abs(theta - link_angle_deg(gw, user)));
        AF = beamforming_fun(link_angle_deg(gw, user)+div_theta(user), N, 0);
        RSS(gw, user) = Tx_power(gw) + AF(idx) + 10 * log10(N) - link_attenuation(gw, user);
    end
end

end

% ------------------------------------------------------------------ %

function [AF_dB] = beamforming_fun(angle_direction_deg, N_array_element, display_mode)
frequency = 915e6; % Frequency in Hz
c = 3e8; % Speed of light in m/s
lambda = c / frequency; % Wavelength in meters
d = lambda / 2; % Distance between antennas (half-wavelength spacing)

% Angles for plotting
theta = linspace(-pi/2, pi/2, 360); % Angles from -90 to 90 degrees in radians

% Steering angle (incoming signal direction)
steering_angle = angle_direction_deg * pi / 180; % 30 degrees converted to radians

% Steering vector
steering_vector = exp(1j * (0:(N_array_element-1))' * (2*pi*d/lambda) * sin(theta - steering_angle));

% Apply Hamming window to the steering vector
window = hanning(N_array_element);
windowed_steering_vector = window .* steering_vector;

% Array factor (Beamforming)
AF = sum(windowed_steering_vector, 1);

% Normalize the array factor
AF_normalized = AF / max(abs(AF));

% Convert to dB scale
AF_dB = 20 * log10(abs(AF_normalized));

% Plot the antenna pattern in polar coordinates
if display_mode == 1 
figure;
polarplot(theta, AF_dB);
rticks([-40 -30 -20 -10 0]); % Define radius ticks in dB
rlim([-40 0]); % Define radius limits
thetalim([-90 90]); % Define theta limits from -90 to 90 degrees
set(gca, 'ThetaZeroLocation', 'top'); % Set zero angle to top
set(gca, 'ThetaDir', 'clockwise'); % Set angle direction to clockwise
grid on;
end
end
