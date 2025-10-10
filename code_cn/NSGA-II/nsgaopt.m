function defaultopt = nsgaopt()
% 函数: defaultopt = nsgaopt()
% 描述: 创建NSGA-II默认选项结构体。
% 语法:  opt = nsgaopt()
%         LSSSSWC, NWPU
%   Revision: 1.3  Data: 2011-07-13
%*************************************************************************


defaultopt = struct(...
... % 优化模型
    'popsize', 50,...           % 种群大小
    'maxGen', 100,...           % 最大代数
    'numVar', 0,...             % 设计变量数量
    'numObj', 0,...             % 目标数量
    'numCons', 0,...            % 约束数量
    'lb', [],...                % 设计变量下界 [1:numVar]
    'ub', [],...                % 设计变量上界 [1:numVar]
    'vartype', [],...           % 变量数据类型 [1:numVar]，1=实数, 2=整数
    'objfun', @objfun,...       % 目标函数
... % 优化模型各组件名称
    'nameObj',{{}},...
    'nameVar',{{}},...
    'nameCons',{{}},...
... % 初始化和输出
    'initfun', {{@initpop}},...         % 种群初始化函数 (默认为随机数)
    'outputfuns',{{@output2file}},...   % 输出函数
    'outputfile', 'populations.txt',... % 输出文件名
    'outputInterval', 1,...             % 输出间隔
    'plotInterval', 5,...               % 两次调用 "plotnsga" 之间的间隔。
... % 遗传算法算子
    'crossover', {{'intermediate', 1.2}},...         % 交叉算子 (比率=1.2)
    'mutation', {{'gaussian',0.1, 0.5}},...          % 变异算子 (缩放因子=0.1, 收缩因子=0.5)
    'crossoverFraction', 'auto', ...                 % 个体中变量的交叉率
    'mutationFraction', 'auto',...                   % 个体中变量的变异率
... % 算法参数
    'useParallel', 'no',...                          % 并行计算种群的目标函数。 {'yes','no'}
    'poolsize', 0,...                                % 并行计算使用的工作进程数，0 = 自动选择。
... % R-NSGA-II 参数
    'refPoints', [],...                              % 用于指定偏好的参考点。每一行是一个参考点。
    'refWeight', [],...                              % 计算欧氏距离时使用的权重因子
    'refUseNormDistance', 'front',...                % 使用由可能的最大和最小目标值归一化的欧氏距离。 {'front','ever','no'}
    'refEpsilon', 0.001 ...                          % 基于epsilon的选择策略中使用的参数
);