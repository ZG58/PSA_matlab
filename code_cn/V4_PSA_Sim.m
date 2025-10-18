% 清理工作区和命令窗口，设置输出格式
clc;
clear;
format long;
rng('shuffle'); % 初始化随机数生成器

% 加载预定义参数
try
    load('Params.mat');
catch
    disp('警告: 未找到 Params.mat 文件。将使用示例数据。');
    IsothermPar = rand(20, 5); % 示例数据
    SimParam = rand(20, 8);    % 示例数据
end

% --- 主要参数设置 ---
N = 10;
type = 'ProcessEvaluation';
simulation_index = 12;

% --- 提取特定材料的参数 ---
IsothermParams     = IsothermPar(simulation_index, :);
material_propertry = SimParam(simulation_index, :);
material    = {};
material{1} = material_propertry;
material{2} = IsothermParams;

% --- 定义随机采样变量 x 的边界 ---
lb = [1e5,  10, 0.01, 0.1, 0, 1e4, 0, 0, 0];      % x 的下界
ub = [10e5, 1000, 0.99, 2, 1, 5e4, 20, 30, 70];     % x 的上界

% lb = [1e5,  10, 0.01, 0.1, 0, 1e4];      % x 的下界
% ub = [10e5, 1000, 0.99, 2, 1, 5e4];     % x 的上界

% --- 在循环开始前，一次性完成所有采样 ---
num_simulations = 1; % 设置您想运行的总仿真次数
%% 

% fprintf('正在为 %d 次仿真预先生成所有随机样本...\n', num_simulations);
% % 生成一个 num_simulations x length(lb) 的矩阵
% % 每一行都是一组独立的随机样本 x
% X_samples = lb + (ub - lb) .* rand(num_simulations, length(lb));
% disp('采样完成。');

%% 
X_samples = [220676, 47.902,	0.236211,	0.16219,	0.529994,	10673.4,	15.4153,	17.1233, 70];

%% 

% --- 预分配空间存储结果 ---
results_objectives = zeros(num_simulations, 2); 
results_constraints = zeros(num_simulations, 3);

% --- 计时开始 ---
tic;

fprintf('\n开始使用预生成的样本进行仿真...\n');
for i = 1:num_simulations
    
    % 从样本矩阵中提取第 i 行作为当前的 x
    x = X_samples(i, :);
    
    fprintf('--> 正在进行第 %d/%d 次仿真...\n', i, num_simulations);

    % 调用仿真函数 (请替换为您的真实函数)
    [objectives, constraints] = V4_PSACycleSimulation(x, material, type, N);
    
    % 存储结果
    results_objectives(i, :) = objectives;
    results_constraints(i, :) = constraints;
end

% --- 计时结束 ---
total_time = toc;

fprintf('\n所有 %d 次仿真已完成！\n', num_simulations);
fprintf('总仿真时间：%.2f 秒\n', total_time);
fprintf('平均每次仿真时间：%.2f 秒\n', total_time / num_simulations);
