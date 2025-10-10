function pop = mutationOp(opt, pop, state)
% 函数: pop = mutationOp(opt, pop, state)
% 描述: 变异算子。所有个体都会进行变异，但
%   一个个体中只有"mutationFraction"比例的设计变量会被改变。
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-13
%*************************************************************************

%*************************************************************************
% 1. 检查参数
%*************************************************************************
% 变异方法
strfun = lower(opt.mutation{1});
numOptions = length(opt.mutation) - 1;
[mutationopt{1:numOptions}] = opt.mutation{2:end};

switch (strfun)
    case 'gaussian'
        fun = @mutationGaussian;
    otherwise
        error('NSGA2:MutationOpError', '不支持的变异算子！');
end

nVar = opt.numVar;

% "auto" 变异率
if( ischar(opt.mutationFraction) )
    if( strcmpi(opt.mutationFraction, 'auto') )
        fraction = 2.0 / nVar;
    else
        error('NSGA2:MutationOpError', '"mutationsFraction" 参数应为标量或 "auto" 字符串。');
    end
else
    fraction = opt.mutationFraction;
end


% 所有个体都会被修改，但一个个体中只有'mutationFraction'比例的设计
% 变量会被改变。
for ind = 1:length(pop)
        child = fun( pop(ind), opt, state, fraction, mutationopt);
        
        % 对整数变量进行四舍五入
        for v = 1:nVar
            if( opt.vartype(v) == 2)
                child.var(v) = round( child.var(v) );
            end
        end

        child.var = varlimit(child.var, opt.lb, opt.ub);
        
        pop(ind) = child;
end



function child = mutationGaussian( parent, opt, state, fraction, options)
% 函数: child = mutationGaussian( parent, opt, state, fraction, options)
% 描述: 高斯变异算子。参考Matlab帮助文档：
%   Genetic Algorithm Options :: Options Reference (Global Optimization Toolbox)
% 参数: 
%   fraction : 个体中变量的变异率
%   options{1} : 缩放因子。对于整数变量，此参数应足够大
%     以使其能从一个值变为另一个值。
%   options{2} : 收缩因子
% 返回: 
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-13
%*************************************************************************


%*************************************************************************
% 1. 验证参数。
%*************************************************************************
if( length(options)~=2)
    error('NSGA2:MutationOpError', '变异算子参数错误！');
end


%*************************************************************************
% 2. 计算"scale"和"shrink"参数。
%*************************************************************************
scale = options{1};
shrink = options{2};
scale = scale - shrink * scale * state.currentGen / opt.maxGen;

lb = opt.lb;
ub = opt.ub;
scale = scale * (ub - lb);


%*************************************************************************
% 3. 执行变异。
%*************************************************************************
child = parent;
numVar = length(child.var);
for i = 1:numVar
    if(rand() < fraction)
        child.var(i) = parent.var(i) + scale(i) * randn();
    end
end