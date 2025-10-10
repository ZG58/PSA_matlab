function opt = verifyOpt(opt)
% 函数: opt = verifyOpt(opt)
% 描述: 验证优化模型。
%         LSSSSWC, NWPU
%   Revision: 1.1  Data: 2011-07-15
%*************************************************************************


%*************************************************************************
%  种群大小 (popsize)
%*************************************************************************
if ( mod(opt.popsize, 2) ~= 0 )
    %warning('NSGA2:PopSizeError', '种群大小应为偶数！%d => %d', opt.popsize, opt.popsize+1);
    opt.popsize = opt.popsize + 1;
end

%*************************************************************************
% 上下界 (lb, ub)
%*************************************************************************
if( length(opt.lb)~=opt.numVar || length(opt.lb)~=opt.numVar )
    error('NSGA2:OptModelError', '上下界的数量(%d,%d)应等于设计变量的数量(%d)！', ...
        length(opt.ub), length(opt.lb), opt.numVar);
end

%*************************************************************************
% 变量类型 (vartype)
%*************************************************************************
if( length(opt.vartype) ~= opt.numVar )
    %warning('NSGA2:OptModelWarning', '设计变量数据类型错误！所有类型都将设为实数编码 (vartype=1)！');
    opt.vartype = ones(1, opt.numVar);
end

%*************************************************************************
% 目标、变量、约束的名称 (nameObj, nameVar, nameCons)
%*************************************************************************
if( ~iscell(opt.nameObj) || ~iscell(opt.nameVar) || ~iscell(opt.nameCons))
    error('NSGA2:OptModelError', '目标、设计变量或约束的名称应在元胞数组中指定，例如 {''obj1'',''obj2''}');
end

if( (~isempty(opt.nameObj)  && length(opt.nameObj)~=opt.numObj) || ...
    (~isempty(opt.nameVar)  && length(opt.nameVar)~=opt.numVar) || ...
    (~isempty(opt.nameCons) && length(opt.nameCons)~=opt.numCons))
    error('NSGA2:OptModelError', '如果指定了任何一个名称，则必须指定所有目标、设计变量或约束的名称！');
end

%*************************************************************************
% 并行计算 (useparallel)
%*************************************************************************
if( ~ischar(opt.useParallel) || ...
    isempty( find(strcmpi(opt.useParallel, {'yes', 'no'}))) )
    error('NSGA2:OptParamError', 'useParallel 只能是 "yes" 或 "no"！');
end

%*************************************************************************
% R-NSGA-II 参数
%*************************************************************************
% 参考点 (refPoints)
if( ~isempty(opt.refPoints) && size(opt.refPoints,2)~=opt.numObj)
    error('NSGA2:OptParamError', '参考点的格式应为 refPoints(nPoint, numObj)！');
end
% 参考点权重 (refWeight)
if( ~isempty(opt.refPoints) && ~isempty(opt.refWeight) && length(opt.refWeight)~=opt.numObj)
    error('NSGA2:OptParamError', 'R-NSGA-II中使用的权重因子向量长度必须等于numObj！');
end