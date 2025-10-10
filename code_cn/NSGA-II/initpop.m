function pop = initpop(opt, pop, varargin)
% 函数: pop = initpop(opt, pop, varargin)
% 描述: 初始化种群。
% 语法:
%   pop = initpop(opt, pop)
%     (默认) 创建一个均匀分布的随机初始种群。
%
%   pop = initpop(opt, pop, 'pop.txt')
%     从已有文件中加载上一次的种群。如果文件中的种群大小
%     小于当前设定的种群大小，则用随机数填充剩余部分。
%
%   pop = initpop(opt, pop, 'pop.txt', ngen)
%     从文件中加载指定代数的种群。
%
%   pop = initpop(opt, pop, oldresult)
%     指定已有的结果结构体。
%
%   pop = initpop(opt, pop, oldresult, ngen)
%     指定已有的结果结构体和将要使用的种群代数。
%
% 参数:
%   pop : 一个空种群
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-01
%*************************************************************************


%*************************************************************************
% 1. 识别参数
%*************************************************************************
method = 'uniform';

if(nargin >= 3)
    if( ischar(varargin{1}) )
        method = 'file';
    elseif( isstruct(varargin{1}) )
        method = 'existpop';
    end
end

%*************************************************************************
% 2. 使用不同方法初始化种群
%*************************************************************************
if( strcmpi(method, 'uniform'))
    pop = initpopUniform(opt, pop);
elseif(strcmpi(method, 'file'))
    fprintf('...从文件 "%s" 初始化种群\n', varargin{1});
    pop = initpopFromFile(opt, pop, varargin{:});
elseif(strcmpi(method, 'existpop'))
    fprintf('...从指定结果初始化种群。\n');
    pop = initpopFromExistResult(opt, pop, varargin{:});
end




function pop = initpopFromFile(opt, pop, varargin)
% 函数: pop = initpopFromFile(opt, pop, varargin)
% 描述: 从指定的种群文件中加载种群。
% 语法:
%   pop = initpop(opt, pop, 'pop.txt')
%   pop = initpop(opt, pop, 'pop.txt', ngen)
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-01
%*************************************************************************
fileName = varargin{1};

oldResult = loadpopfile(fileName);
pop = initpopFromExistResult(opt, pop, oldResult, varargin{2:end});





function pop = initpopFromExistResult(opt, pop, varargin)
% 函数: pop = initpopFromExistResult(opt, pop, varargin)
% 描述: 从已有的结果结构体中加载种群。
% 语法:
%   pop = initpop(opt, pop, oldresult)
%   pop = initpop(opt, pop, oldresult, ngen)
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-01
%*************************************************************************

% 1. 验证参数
oldresult = varargin{1};
if( ~isstruct(oldresult) || ~isfield(oldresult, 'pops') )
    error('NSGA2:InitPopError', '指定的结果结构体不正确！');
end


oldpops = oldresult.pops;
ind = oldpops(1,1);     % 用于验证优化参数的个体
if( opt.numVar ~= length(ind.var) || ...
    opt.numObj ~= length(ind.obj) || ...
    opt.numCons ~= length(ind.cons) )
    error('NSGA2:InitPopError', ...
        '指定的优化结果与当前优化模型不匹配！');
end
clear ind


% 2. 决定使用哪一代种群
ngen = 0;
if( nargin >= 4)
    ngen = varargin{2};
end

maxGen = size(oldpops, 1);
if(ngen == 0)
    ngen = maxGen;
elseif(ngen > maxGen)
    warning('NSGA2:InitPopWarning', ...
        '指定的代数 "%d" 不存在, 使用 "%d" 代替。',...
        ngen, maxGen);
    ngen = maxGen;
end


% 3. 创建初始种群
popsizeOld = size(oldpops, 2);
popsizeNew = opt.popsize;

if( popsizeNew <= popsizeOld )      % a) 全部来自旧种群
    for i = 1:popsizeNew
        pop(i).var = oldpops(ngen, i).var;
    end
else                                % b) 使用随机个体填充种群
    for i = 1:popsizeOld
        pop(i).var = oldpops(ngen, i).var;
    end
    pop(popsizeOld+1:end) = initpopUniform(opt, pop(popsizeOld+1:end));
end




function pop = initpopUniform(opt, pop)
% 函数: pop = initpopUniform(opt, pop)
% 描述: 使用随机数初始化种群
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-01
%*************************************************************************

nVar = opt.numVar;
type = opt.vartype;

lb = opt.lb;
ub = opt.ub;

popsize = length(pop);
for i = 1:popsize
    var = lb + rand(1, nVar) .* (ub-lb);
    
    % 如果设计变量是整数，则四舍五入到最近的整数
    for v = 1:nVar
        if( type(v) == 2)
            var(v) = round(var(v));
        end
    end
    
    % 限制在上下界内
    var = varlimit(var, lb, ub);

    pop(i).var = var;

end