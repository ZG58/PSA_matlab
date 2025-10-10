clc
clear
format long
delete(gcp('nocreate'));
parpool('local', 4);

addpath('CycleSteps')
addpath('NSGA-II')

load('Params')

N = 10 ;
type = 'EconomicEvaluation' ;

for i = 12:12
    
% 加载参数
IsothermParams     = IsothermPar(i, :) ;
material_propertry = SimParam(i, :)    ;

material    = {}                 ;
material{1} = material_propertry ;
material{2} = IsothermParams     ;

Function = @(x) PSACycleSimulation( x, material, type, N ) ; % 用于模拟PSA循环的函数

% 初始变量
[~, vars] = sortt(loadpopfile('UTSA-16_Process.txt'));
% vars = [vars, ones(length(vars), 1), 1e4*ones(length(vars), 1)];

options            = nsgaopt();                          % 创建默认选项结构体
options.popsize    = 60;                                 % 种群大小
options.outputfile = 'UTSA-16_Economic.txt';
options.maxGen     = 80;                                 % 最大代数

options.vartype    = [1, 1, 1, 1, 1, 1] ;

options.initfun={@Pop_Override, vars}   ;                % 从先前结果提供变量

options.numObj  = 2 ;                                    % 目标数量
options.numVar  = 6 ;                                    % 设计变量数量
options.numCons = 3 ;                                    % 约束数量
options.lb = [1e5,  10, 0.01, 0.1, 0, 1e4]   ;           % x的下界
options.ub = [10e5, 1000, 0.99, 2, 1, 5e4]   ;           % x的上界
options.nameObj = {'-productivity','energy'} ;           % 目标名称会显示在GUI窗口中。
options.objfun  = Function                   ;           % 目标函数句柄

options.useParallel = 'yes' ;                            % 此处并行计算不是必需的
options.poolsize     =  4   ;                            % 工作进程数

result = nsga2(options)     ;                            % 开始优化！


% 使用30个有限体积单元重新优化
N = 30 ; 
Function = @(x) FiveStepModSkarstromProcessSim( x, material, type, N ) ; % 用于模拟PSA循环的函数

options.objfun     = Function           ;                % 目标函数句柄
options.initfun    = {@initpop, result} ;                % 从先前结果提供变量
options.maxGen     = 120                ;                % 种群大小
options.outputfile = 'UTSA-16_Economic_2.txt'         ;
result2            = nsga2(options)     ;                % 开始优化！

end