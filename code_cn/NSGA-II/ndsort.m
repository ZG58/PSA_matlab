function [opt, pop] = ndsort(opt, pop)
% 函数: [opt, pop] = ndsort(opt, pop)
% 描述: 快速非支配排序。
%
%         LSSSSWC, NWPU
%    Revision: 1.4  Data: 2011-07-26
%*************************************************************************




%*************************************************************************
% 1. 初始化变量
%   indi.np：支配该个体的个体数量
%   indi.sp(:): 该个体所支配的个体集合
%*************************************************************************
N = length(pop);    %种群大小
ind = repmat(struct('np',0, 'sp', []),[1,N]);

for i = 1:N
    pop(i).rank = 0;
    pop(i).distance = 0;
    pop(i).prefDistance = 0;
end


%*************************************************************************
% 2. 快速非支配排序
%*************************************************************************
% 计算支配矩阵以提高效率。

% 注意: 在我的电脑上，Matlab 2010b版本中，"for"语句比"vertcat"语句效率更高。我不知道为什么。(LSSSSWC, 2011-07-25)
nViol   = zeros(N, 1);
violSum = zeros(N, 1);
for i = 1:N
    nViol(i)    = pop(i).nViol;
    violSum(i)  = pop(i).violSum;
end
% nViol   = vertcat(pop(:).nViol);
% violSum = vertcat(pop(:).violSum);

obj     = vertcat(pop(:).obj);
domMat  = calcDominationMatrix(nViol, violSum, obj); % 用于提高效率的支配矩阵


% 计算每个个体的np和sp
for p = 1:N-1
    for q = p+1:N
        if(domMat(p, q) == 1)          % p 支配 q
            ind(q).np = ind(q).np + 1;
            ind(p).sp = [ind(p).sp , q];
        elseif(domMat(p, q) == -1)     % q 支配 p
            ind(p).np = ind(p).np + 1;
            ind(q).sp = [ind(q).sp , p];
        end
    end
end


% 第一个前沿（等级 = 1）
front(1).f = [];    % 'front'结构体中只有一个字段'f'。
                    % 这是故意的，因为每个前沿中的
                    % 个体数量是不同的。
for i = 1:N
    if( ind(i).np == 0 )
        pop(i).rank = 1;
        front(1).f = [front(1).f, i];
    end
end

% 计算每个个体的帕累托等级，即 pop(:).rank
fid = 1;        % 帕累托前沿ID
while( ~isempty(front(fid).f) )
    Q = [];
    for p = front(fid).f
        for q = ind(p).sp
            ind(q).np = ind(q).np -1;
            if( ind(q).np == 0 )
                pop(q).rank = fid+1;
                Q = [Q, q];
            end
        end
    end
    fid = fid + 1;
    
    front(fid).f = Q;
end
front(fid) = [];    % 删除最后一个空的前沿集合



%*************************************************************************
% 3. 计算距离
%*************************************************************************
if(isempty(opt.refPoints))
    pop = calcCrowdingDistance(opt, pop, front);
else
    [opt, pop] = calcPreferenceDistance(opt, pop, front);
end





function domMat = calcDominationMatrix(nViol, violSum, obj)
% 函数: domMat = calcDominationMatrix(nViol, violSum, obj)
% 描述: 计算支配矩阵，该矩阵使用约束支配关系指定两个个体
%   之间的支配关系。
%
% 返回: 
%   domMat(N,N) : 支配矩阵
%       domMat(p,q)=1  : p 支配 q
%       domMat(p,q)=-1 : q 支配 p
%       domMat(p,q)=0  : 互不支配
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-13
%*************************************************************************

N       = size(obj, 1);
numObj  = size(obj, 2);

domMat  = zeros(N, N);

for p = 1:N-1
    for q = p+1:N
        %*************************************************************************
        % 1. p 和 q 都是可行的
        %*************************************************************************
        if(nViol(p) == 0 && nViol(q)==0)
            pdomq = false;
            qdomp = false;
            for i = 1:numObj
                if( obj(p, i) < obj(q, i) )         % 目标函数是最小化！
                    pdomq = true;
                elseif(obj(p, i) > obj(q, i))
                    qdomp = true;
                end
            end

            if( pdomq && ~qdomp )
                domMat(p, q) = 1;
            elseif(~pdomq && qdomp )
                domMat(p, q) = -1;
            end
        %*************************************************************************
        % 2. p 可行, q 不可行
        %*************************************************************************
        elseif(nViol(p) == 0 && nViol(q)~=0)
            domMat(p, q) = 1;
        %*************************************************************************
        % 3. q 可行, p 不可行
        %*************************************************************************
        elseif(nViol(p) ~= 0 && nViol(q)==0)
            domMat(p, q) = -1;
        %*************************************************************************
        % 4. p 和 q 都是不可行的
        %*************************************************************************
        else
            if(violSum(p) < violSum(q))
                domMat(p, q) = 1;
            elseif(violSum(p) > violSum(q))
                domMat(p, q) = -1;
            end
        end
    end
end

domMat = domMat - domMat';





function [opt, pop] = calcPreferenceDistance(opt, pop, front)
% 函数: [opt, pop] = calcPreferenceDistance(opt, pop, front)
% 描述: 计算R-NSGA-II中使用的'偏好距离'。
% 返回: 
%   opt : 仅当 opt.refUseNormDistance=='ever' 时，此结构体可能被修改。
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.1  Data: 2011-07-26
%*************************************************************************

%*************************************************************************
% 1. 初始化
%*************************************************************************
numObj = length( pop(1).obj );  % 目标数量

refPoints = opt.refPoints;
refWeight = opt.refWeight;      % 目标权重因子
if(isempty(refWeight))
    refWeight = ones(1, numObj);
end
epsilon = opt.refEpsilon;
numRefPoint = size(refPoints, 1);

% 决定归一化因子
bUseFrontMaxMin = false;    % bUseFrontMaxMin : 是否使用前沿中的最大和最小值作为归一化因子。
if( strcmpi(opt.refUseNormDistance, 'ever') )
    % 1) 找到每个目标可能的最大和最小值（不限于当前种群）
    obj = vertcat(pop.obj);
    if( ~isfield(opt, 'refObjMax_tmp') )
        opt.refObjMax_tmp = max(obj);
        opt.refObjMin_tmp = min(obj);
    else
        objMax = max(obj);
        objMin = min(obj);
        for i = 1:numObj
            if(opt.refObjMax_tmp(i) < objMax(i))
                opt.refObjMax_tmp(i) = objMax(i);
            end
            if(opt.refObjMin_tmp(i) > objMin(i))
                opt.refObjMin_tmp(i) = objMin(i);
            end
        end
        clear objMax objMin
    end
    objMaxMin = opt.refObjMax_tmp - opt.refObjMin_tmp;
    clear obj
elseif( strcmpi(opt.refUseNormDistance, 'front') )
    % 2) 不使用归一化的欧氏距离。
    bUseFrontMaxMin = true;
elseif( strcmpi(opt.refUseNormDistance, 'no') )
    % 3) 不使用归一化的欧氏距离。
    objMaxMin = ones(1,numObj);
else
    % 3) 错误
    error('NSGA2:ParamError', ...
        '不支持的参数: options.refUseNormDistance="%s", 只支持 "yes" 或 "no"',...
        opt.refUseNormDistance);
end


%*************************************************************************
% 2. 计算偏好距离 pop(:).prefDistance
%*************************************************************************
for fid = 1:length(front)
    % 步骤1: 计算每个前沿中的加权欧氏距离
    idxFront = front(fid).f;            % idxFront : 当前前沿中个体的索引
    numInd = length(idxFront);          % numInd : 当前前沿中的个体数量
    popFront = pop(idxFront);           % popFront : 在前沿fid中的个体

    objFront = vertcat(popFront.obj);   % objFront : 所有个体的全部目标值

    if(bUseFrontMaxMin)
        objMaxMin = max(objFront) - min(objFront); % objMaxMin : 当前前沿中的归一化因子
    end

    % normDistance : 加权归一化欧氏距离
    normDistance = calcWeightNormDistance(objFront, refPoints, objMaxMin, refWeight);
    
    
    % 步骤2: 分配偏好距离
    prefDistanceMat = zeros(numInd, numRefPoint);
    for ipt = 1:numRefPoint
        [~,ix] = sort(normDistance(:, ipt));
        prefDistanceMat(ix, ipt) = 1:numInd;
    end
    prefDistance = min(prefDistanceMat, [], 2);
    clear ix

    
    % 步骤3: Epsilon 清除策略
    idxRemain = 1:numInd;           % idxRemain : 未处理个体的索引
    while(~isempty(idxRemain))
        % 1. 从剩余个体中选择一个
        objRemain = objFront( idxRemain, :);
        selIdx = randi( [1,length(idxRemain)] );
        selObj = objRemain(selIdx, :);

        % 2. 计算归一化欧氏距离
        % distanceToSel : 到所选点的归一化欧氏距离
        distanceToSel = calcWeightNormDistance(objRemain, selObj, objMaxMin, refWeight);
        

        % 3. 处理epsilon邻域内的个体
        idx = find( distanceToSel <= epsilon );     % idx : 在 idxRemain 中的索引
        if(length(idx) == 1)    % 唯一的个体就是被选中的那个
            idxRemain(selIdx)=[];
        else
            for i=1:length(idx)
                if( idx(i)~=selIdx )
                    idInIdxRemain = idx(i);     % idx 是在 idxRemain 向量中的索引
                    id = idxRemain(idInIdxRemain);
                    
                    % *增加偏好距离以降低这些个体
                    % 被选择的几率。
                    prefDistance(id) = prefDistance(id) + round(numInd/2);
                end
            end
            idxRemain(idx) = [];
        end
        
    end

    % 保存偏好距离
    for i=1:numInd
        id = idxFront(i);
        pop(id).prefDistance = prefDistance(i);
    end
end


function distance = calcWeightNormDistance(points, refPoints, maxMin, weight)
% 函数: calcWeightNormDistance(points, refPoints, maxMin, weight)
% 描述: 计算从"points"到"refPoints"的加权欧氏距离
% 参数: 
%   points(nPoint, N)       : 每一行是N维空间中的一个点。
%   refPoints(nRefPoint, N) : 每一行是一个参考点。
%   maxMin(1, N)            : 归一化因子。
%   weight(1, N)            : 权重
%
% 返回: 
%   distance(nPoint, nRefPoint)
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-14
%*************************************************************************

nRefPoint = size(refPoints, 1);     % 参考点数量
nPoint = size(points, 1);           % 点的数量

distance = zeros(nPoint, nRefPoint);
for ipt = 1:nRefPoint
    refpt = refPoints(ipt, :);
    for i = 1:nPoint
        weightNormDist = ((points(i, :)-refpt) ./ maxMin).^2 .* weight;
        distance(i, ipt) = sqrt(sum(weightNormDist));
    end
end





function pop = calcCrowdingDistance(opt, pop, front)
% 函数: pop = calcCrowdingDistance(opt, pop, front)
% 描述: 计算原始NSGA-II中使用的'拥挤度距离'。
% 语法:
% 参数: 
% 返回: 
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-11
%*************************************************************************

numObj = length( pop(1).obj );  % 目标数量
for fid = 1:length(front)
    idx = front(fid).f;
    frontPop = pop(idx);        % frontPop : 在前沿fid中的个体
    
    numInd = length(idx);       % nInd : 当前前沿中的个体数量
    
    obj = vertcat(frontPop.obj);
    obj = [obj, idx'];          % 目标值与个体ID一起排序
    for m = 1:numObj
        obj = sortrows(obj, m);

        colIdx = numObj+1;
        pop( obj(1, colIdx) ).distance = Inf;         % 第一个个体
        pop( obj(numInd, colIdx) ).distance = Inf;    % 最后一个个体
        
        minobj = obj(1, m);         % 目标m的最小值
        maxobj = obj(numInd, m);    % 目标m的最大值
        
        for i = 2:(numInd-1)
            id = obj(i, colIdx);
            pop(id).distance = pop(id).distance + (obj(i+1, m) - obj(i-1, m)) / (maxobj - minobj);
        end
    end
end