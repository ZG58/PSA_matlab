function pop= Pop_Override(opt, pop, varargin)
% 函数: pop= Pop_Override(opt, pop, varargin)
% 描述: 从已有变量中加载种群。不会因先前结果的
% 约束数量不匹配而拒绝

    popsize_new=opt.popsize;

    variable=opt.initfun{2};
    
    popsize_old=size(variable, 1);

    numb_var_sim=opt.numVar;

    numb_var_supp=size(variable, 2);

    if numb_var_sim ~= numb_var_supp
        error('NSGA2:OptModelError', '提供的变量数量与要求的数量不相等');
    end

    pop_size=opt.popsize;

    variable=opt.initfun{2};
    
    if popsize_old >= popsize_new
        for j = 1:pop_size
            pop(j).var = variable(j, :);
        end
    else
        for j = 1:popsize_old
            pop(j).var = variable(j, :);
        end
        pop(popsize_old+1:end) = initpopUniform(opt, pop(popsize_old+1:end));
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
    end

end