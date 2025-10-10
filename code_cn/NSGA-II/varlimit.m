function var = varlimit(var, lb, ub)
% 函数: var = varlimit(var, lb, ub)
% 描述: 将变量限制在[lb, ub]范围内。
%
%         LSSSSWC, NWPU
%    Revision: 1.0  Data: 2011-04-20
%*************************************************************************

numVar = length(var);
for i = 1:numVar
    if( var(i) < lb(i) )
        var(i) = lb(i);
    elseif( var(i) > ub(i) )
        var(i) = ub(i);
    end
end