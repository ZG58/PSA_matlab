function state = statpop(pop, state)
% 函数: state = statpop(pop, state)
% 描述: 统计种群信息。
%
%         LSSSSWC, NWPU
%    Revision: 1.0  Data: 2011-04-20
%*************************************************************************


N = length(pop);
rankVec = vertcat(pop.rank);
rankVec = sort(rankVec);

state.frontCount = rankVec(N);
state.firstFrontCount = length( find(rankVec==1) );