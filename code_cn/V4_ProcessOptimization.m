clc
format long
%delete(gcp('nocreate'));
%parpool('local', 4);

addpath('CycleSteps')
addpath('NSGA-II')

load('Params')

N = 10 ;
type = 'ProcessEvaluation' ;

for i = 12:12
    
% 加载参数
IsothermParams     = IsothermPar(i, :) ;
material_propertry = SimParam(i, :)    ;

material    = {}                 ;
material{1} = material_propertry ;
material{2} = IsothermParams     ;

Function = @(x) V4_PSACycleSimulation( x, material, type, N ) ; % 用于模拟PSA循环的函数

options         = nsgaopt() ;                            % 创建默认选项结构体
options.popsize = 8        ;                            % 种群大小
options.maxGen  = 1        ;                            % 最大代数

options.vartype    = [1, 1, 1, 1, 1, 1, 1, 1]         ;
options.outputfile = 'UTSA-16_Process.txt' ;

options.numObj  = 2 ;                                    % 目标数量
options.numVar  = 8 ;                                    % 设计变量数量
options.numCons = 3 ;                                    % 约束数量

options.lb = [1e5,  10, 0.01, 0.1, 0, 1e4, 0, 0];                % x 的下界
options.ub = [10e5, 1000, 0.99, 2, 1, 5e4, 20, 30];              % x 的上界

options.nameObj = {'-purity','recovery'} ;               % 目标名称会显示在GUI窗口中。
options.objfun  = Function               ;               % 目标函数句柄

options.useParallel = 'no' ;                             % 此处并行计算不是必需的
%options.poolsize     = 4   ;                            % 工作进程数

result = nsga2(options)     ;                            % 开始优化！


end