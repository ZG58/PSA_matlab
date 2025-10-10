function newpop = selectOp(opt, pop)
% 函数: newpop = selectOp(opt, pop)
% 描述: 选择算子，使用二元锦标赛选择法。
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-12
%*************************************************************************

popsize = length(pop);
pool = zeros(1, popsize);   % pool : 被选中的个体索引

randnum = randi(popsize, [1, 2 * popsize]);

j = 1;
for i = 1:2:(2*popsize)
    p1 = randnum(i);
    p2 = randnum(i+1);
    
    if(~isempty(opt.refPoints))
        % 偏好算子 (R-NSGA-II)
        result = preferenceComp( pop(p1), pop(p2) );
    else
        % 拥挤度比较算子 (NSGA-II)
        result = crowdingComp( pop(p1), pop(p2) );
    end
    
    if(result == 1)
        pool(j) = p1;
    else
        pool(j) = p2;
    end
    
    j = j + 1;
end
newpop = pop( pool );



function result = crowdingComp( guy1, guy2)
% 函数: result = crowdingComp( guy1, guy2)
% 描述: 拥挤度比较算子。
% 返回: 
%   1 = guy1 优于 guy2
%   0 = 其他情况
%
%         LSSSSWC, NWPU
%    Revision: 1.0  Data: 2011-04-20
%*************************************************************************

if((guy1.rank < guy2.rank) || ((guy1.rank == guy2.rank) && (guy1.distance > guy2.distance) ))
    result = 1;
else
    result = 0;
end



function result = preferenceComp(guy1, guy2)
% 函数: result = preferenceComp(guy1, guy2)
% 描述: R-NSGA-II中使用的偏好算子
% 返回: 
%   1 = guy1 优于 guy2
%   0 = 其他情况
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-11
%*************************************************************************

if(  (guy1.rank  < guy2.rank) || ...
    ((guy1.rank == guy2.rank) && (guy1.prefDistance < guy2.prefDistance)) )
    result = 1;
else
    result = 0;
end