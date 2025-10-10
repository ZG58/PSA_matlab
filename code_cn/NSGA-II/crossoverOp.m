function pop = crossoverOp(opt, pop, state)
% 函数: pop = crossoverOp(opt, pop, state)
% 描述: 交叉算子。所有个体都会进行交叉操作，但是
%   一个个体中只有"crossoverFraction"比例的设计变量会被改变。
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-13
%*************************************************************************

%*************************************************************************
% 1. 检查参数
%*************************************************************************
% 决定交叉方法
strfun = lower(opt.crossover{1});
numOptions = length(opt.crossover) - 1;
[crossoverOpt{1:numOptions}] = opt.crossover{2:end};

switch( strfun )
    case 'intermediate'
        fun = @crsIntermediate;
    otherwise
        error('NSGA2:CrossoverOpError', '不支持的交叉算子！');
end

nVar = opt.numVar;

% "auto" 交叉率
if( ischar(opt.crossoverFraction) )
    if( strcmpi(opt.crossoverFraction, 'auto') )
        fraction = 2.0 / nVar;
    else
        error('NSGA2:CrossoverOpError', '"crossoverFraction" 参数应为标量或 "auto" 字符串。');
    end
else
    fraction = opt.crossoverFraction;
end


for ind = 1:2:length(pop)    % 种群大小应为偶数
        % 创建子代
        [child1, child2] = fun( pop(ind), pop(ind+1), fraction, crossoverOpt );
        
        % 取整
        for v = 1:nVar
            if( opt.vartype(v) == 2)
                child1.var(v) = round( child1.var(v) );
                child2.var(v) = round( child2.var(v) );
            end
        end

        % 边界限制
        child1.var = varlimit(child1.var, opt.lb, opt.ub);
        child2.var = varlimit(child2.var, opt.lb, opt.ub);
        
        pop(ind)     = child1;
        pop(ind+1)   = child2;
    
end



function [child1, child2] = crsIntermediate(parent1, parent2, fraction, options)
% 函数: [child1, child2] = crsIntermediate(parent1, parent2, fraction, options)
% 描述: (用于实数编码) 中间交叉。（与Matlab的交叉算子相同）
%       child = parent1 + rand * Ratio * ( parent2 - parent1)
% 参数: 
%   fraction : 个体中变量的交叉率
%   options = 比率
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-13
%*************************************************************************


if( length(options)~=1 || ~isnumeric(options{1}))
    error('NSGA2:CrossoverOpError', '交叉算子参数错误！');
end

ratio = options{1};

child1 = parent1;
child2 = parent2;

nVar = length(parent1.var);
crsFlag = rand(1, nVar) < fraction;

randNum = rand(1,nVar);     % 均匀分布

child1.var = parent1.var + crsFlag .* randNum .* ratio .* (parent2.var - parent1.var);
child2.var = parent2.var - crsFlag .* randNum .* ratio .* (parent2.var - parent1.var);