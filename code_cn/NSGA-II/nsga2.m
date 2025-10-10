function result = nsga2(opt, varargin)
% 函数: result = nsga2(opt, varargin)
% 描述: NSGA-II的主流程。注意：
%   所有目标函数都必须是最小化问题。如果某个目标是最大化问题，
%   则应将该目标乘以-1。
%
% 语法:
%   result = nsga2(opt): 'opt' 是由函数 nsgaopt() 生成的。
%   result = nsga2(opt, param): 'param' 可以是任何数据类型，它将被
%       传递给目标函数 objfun()。
%
%   然后，结果结构体可以传递给 plotnsga 来显示
%   种群： plotnsga(result);
%
% 参数:
%   opt : 由函数 nsgaopt() 生成的结构体。
%   varargin : 将传递给目标函数的附加参数。
%       它可以是任何数据类型。例如，如果你调用：nsga2(opt, param)，
%       那么 objfun 将被调用为 objfun(x,param)，其中 x 是
%       设计变量向量。
% 返回:
%   result : 包含优化结果的结构体。
%
%         LSSSSWC, NWPU
%   Revision: 1.2  Data: 2011-07-26
%*************************************************************************


tStart = tic();
%*************************************************************************
% 验证优化模型
%*************************************************************************
opt = verifyOpt(opt);

%*************************************************************************
% 变量初始化
%*************************************************************************
nVar    = opt.numVar;
nObj    = opt.numObj;
nCons   = opt.numCons;
popsize = opt.popsize;

% pop : 当前种群
% newpop : 由遗传算子创建的新种群
% combinepop = pop + newpop;
pop = repmat( struct(...
    'var', zeros(1,nVar), ...
    'obj', zeros(1,nObj), ...
    'cons', zeros(1,nCons),...
    'rank', 0,...
    'distance', 0,...
    'prefDistance', 0,...       % R-NSGA-II中使用的偏好距离
    'nViol', 0,...
    'violSum', 0),...
    [1,popsize]);

% state: 某一代替换状态
state = struct(...
'currentGen', 1,...         % 当前代数
'evaluateCount', 0,...      % 目标函数评估次数
'totalTime', 0,...          % 从开始至今的总时间
'firstFrontCount', 0,...    % 第一个前沿的个体数量
'frontCount', 0,...         % 前沿数量
'avgEvalTime', 0 ...        % 目标函数的平均评估时间（当前代）
);

result.pops     = repmat(pop, [opt.maxGen, 1]);     % 每一行是某一代的种群
result.states   = repmat(state, [opt.maxGen, 1]);   % 每一行是某一代的优化状态
result.opt      = opt;                              % 用于输出

% 全局变量
global STOP_NSGA;   % STOP_NSGA : 用于GUI，如果STOP_NSGA~=0，则停止优化
STOP_NSGA = 0;


%*************************************************************************
% 初始化P0种群
%*************************************************************************
ngen = 1;
pop = opt.initfun{1}(opt, pop, opt.initfun{2:end});
[pop, state] = evaluate(opt, pop, state, varargin{:});
[opt, pop] = ndsort(opt, pop);

% 状态
state.currentGen = ngen;
state.totalTime = toc(tStart);
state = statpop(pop, state);

result.pops(1, :) = pop;
result.states(1)  = state;

% 输出
opt = callOutputfuns(opt, state, pop);


%*************************************************************************
% NSGA2迭代
%*************************************************************************
while( ngen < opt.maxGen && STOP_NSGA==0)
    % 0. 显示一些信息
	ngen = ngen+1;
    state.currentGen = ngen;
    
    % 1. 创建新种群
    newpop = selectOp(opt, pop);
    newpop = crossoverOp(opt, newpop, state);
    newpop = mutationOp(opt, newpop, state);
    [newpop, state] = evaluate(opt, newpop, state, varargin{:});

    % 2. 合并新旧种群 : combinepop = pop + newpop
    combinepop = [pop, newpop];
    
    % 3. 快速非支配排序
    [opt, combinepop] = ndsort(opt, combinepop);
    
    % 4. 提取下一代种群
    pop = extractPop(opt, combinepop);

    % 5. 保存当前代的结果
    state.totalTime = toc(tStart);
    state = statpop(pop, state);
    
    result.pops(ngen, :) = pop;
    result.states(ngen)  = state;

    % 6. 输出
    if( mod(ngen, opt.outputInterval)==0 )
        opt = callOutputfuns(opt, state, pop);
    end
    
end

% 调用输出函数以关闭文件
opt = callOutputfuns(opt, state, pop, -1);



%toc(tStart);