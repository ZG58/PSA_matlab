clc
format long
delete(gcp('nocreate'));
parpool('local', 3);

addpath('CycleSteps')
addpath('NSGA-II')

load('Params_2.mat')

N = 10 ;
type = 'ProcessEvaluation' ;

% --- 创建用于存放结果的文件夹 ---
resultsDir = 'Results';   % 定义结果文件夹的名称
if ~exist(resultsDir, 'dir') % 检查名为'Results'的文件夹是否存在
   mkdir(resultsDir);      % 如果不存在，则创建该文件夹
end
% ------------------------------------

for i = 16:17
    
% 加载参数
IsothermParams     = IsothermPar(i, :) ;
material_propertry = SimParam(i, :)    ;

material    = {}                 ;
material{1} = material_propertry ;
material{2} = IsothermParams     ;

Function = @(x) V3_PSACycleSimulation( x, material, type, N ) ; % 用于模拟PSA循环的函数

options         = nsgaopt() ;                            % 创建默认选项结构体
options.popsize = 9        ;                            % 种群大小
options.maxGen  = 60        ;                            % 最大代数

options.vartype    = [1, 1, 1, 1, 1, 1, 1, 1, 1]         ;

% 为每一次循环创建一个唯一的文件名，例如 'Process_Output_1.txt'
outputFileName = sprintf('Process_Output_%d.txt', i);
% 使用fullfile函数构建完整的文件路径，确保跨操作系统兼容性
options.outputfile = fullfile(resultsDir, outputFileName);

options.numObj  = 2 ;                                    % 目标数量
options.numVar  = 9 ;                                    % 设计变量数量
options.numCons = 3 ;                                    % 约束数量

options.lb = [1e5,  10, 0.01, 0.1, 0, 1e4, 0, 0, 0];     % x 的下界
options.ub = [10e5, 1000, 0.99, 2, 1, 5e4, 20, 30, 70];  % x 的上界

options.nameObj = {'-purity','recovery'} ;               % 目标名称会显示在GUI窗口中。
options.objfun  = Function               ;               % 目标函数句柄

options.useParallel = 'yes' ;                             % 此处并行计算不是必需的
options.poolsize     = 3   ;                                % 工作进程数

result = nsga2(options)     ;                            % 开始优化！


end