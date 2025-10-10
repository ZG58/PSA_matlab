function [pop, state] = evaluate(opt, pop, state, varargin)
% 函数: [pop, state] = evaluate(opt, pop, state, varargin)
% 描述: 评估种群中每个个体的目标函数值。
%
%         LSSSSWC, NWPU
%    Revision: 1.0  Data: 2011-04-20
%*************************************************************************

N = length(pop);
allTime = zeros(N, 1);  % allTime : 用于计算平均评估时间

%*************************************************************************
% 并行评估目标函数
%*************************************************************************
if( strcmpi(opt.useParallel, 'yes') == 1 )
%     curPoolInfo = gcp('nocreate');
%     
%     % 检查并行计算池（parpool）是否已开启
%     [Pooldatasize, ~]=size(curPoolInfo);
%     if Pooldatasize == 0
%         curPoolsize=0;
%     else
%         curPoolsize = curPoolInfo.NumWorkers;
%     end
% 
%     % 当前没有开启的工作进程
%     if(curPoolsize == 0)
%         if(opt.poolsize == 0)
%             parpool open local;
%         else
%             parpool(opt.poolsize);
%         end
%     % 关闭并重建工作进程
%     else
%         if(opt.poolsize ~= curPoolsize)
%             delete(gcp);
%             parpool(opt.poolsize);
%         end
%     end

    parfor i = 1:N
        fprintf('\n正在评估目标函数... 第 %d 代 / 共 %d 代, 个体 %d / 共 %d \n',state.currentGen, opt.maxGen, i, N);
        [pop(i), allTime(i)] = evalIndividual(pop(i), opt.objfun, varargin{:});
    end

%*************************************************************************
% 串行评估目标函数
%*************************************************************************
else
    for i = 1:N
        fprintf('\n正在评估目标函数... 第 %d 代 / 共 %d 代, 个体 %d / 共 %d \n',state.currentGen, opt.maxGen, i, N);
        [pop(i), allTime(i)] = evalIndividual(pop(i), opt.objfun, varargin{:});
    end
end

%*************************************************************************
% 统计
%*************************************************************************
state.avgEvalTime   = sum(allTime) / length(allTime);
state.evaluateCount = state.evaluateCount + length(pop);




function [indi, evalTime] = evalIndividual(indi, objfun, varargin)
% 函数: [indi, evalTime] = evalIndividual(indi, objfun, varargin)
% 描述: 评估单个目标函数。
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-25
%*************************************************************************

tStart = tic;
[y, cons] = objfun( indi.var, varargin{:} );
evalTime = toc(tStart);

% 保存目标函数值和约束违反情况
indi.obj = y;

if( ~isempty(indi.cons) )
    idx = find( cons );
    indi.cons=cons;
    if( ~isempty(idx) )
        indi.nViol = length(idx);
        indi.violSum = sum( abs(cons) );
    else
        indi.nViol = 0;
        indi.violSum = 0;
    end

end