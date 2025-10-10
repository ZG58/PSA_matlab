function nextpop = extractPop(opt, combinepop)
% 函数: nextpop = extractPop(opt, combinepop)
% 描述: 从'combinepop'（种群大小为2n）中提取最优的 n 个个体。
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-12
%*************************************************************************

popsize = length(combinepop) / 2;
nextpop = combinepop(1:popsize);    % 仅用于初始化

rankVector = vertcat(combinepop.rank);

n = 0;          % 下一代种群的个体数量
rank = 1;       % 当前等级
idx = find(rankVector == rank);
numInd = length(idx);       % 当前前沿中的个体数量
while( n + numInd <= popsize )
    nextpop( n+1 : n+numInd ) = combinepop( idx );
    
    n = n + numInd;
    rank = rank + 1;
    
    idx = find(rankVector == rank);
    numInd = length(idx);
end

% 如果当前前沿的个体数量加上下一前沿的个体数量
% 大于种群大小，则通过拥挤度距离（NSGA-II）或偏好距离（R-NSGA-II）
% 来选择最优个体。
if( n < popsize )
    if(~isempty(opt.refPoints))
        prefDistance = vertcat(combinepop(idx).prefDistance);
        prefDistance = [prefDistance, idx];
        prefDistance = sortrows( prefDistance, 1);
        idxSelect  = prefDistance( 1:popsize-n, 2);       % 选择具有最小偏好距离的个体
        nextpop(n+1 : popsize) = combinepop(idxSelect);
    else
        distance = vertcat(combinepop(idx).distance);
        distance = [distance, idx];
        distance = flipud( sortrows( distance, 1) );      % 在前沿中按拥挤度距离的降序对个体进行排序。
        idxSelect  = distance( 1:popsize-n, 2);           % 选择（popsize-n）个具有最大拥挤度距离的个体。
        nextpop(n+1 : popsize) = combinepop(idxSelect);
    end
end