function [opt, pop] = ndsort(opt, pop)
% ����: [opt, pop] = ndsort(opt, pop)
% ����: ���ٷ�֧������
%
%         LSSSSWC, NWPU
%    Revision: 1.4  Data: 2011-07-26
%*************************************************************************




%*************************************************************************
% 1. ��ʼ������
%   indi.np��֧��ø���ĸ�������
%   indi.sp(:): �ø�����֧��ĸ��弯��
%*************************************************************************
N = length(pop);    %��Ⱥ��С
ind = repmat(struct('np',0, 'sp', []),[1,N]);

for i = 1:N
    pop(i).rank = 0;
    pop(i).distance = 0;
    pop(i).prefDistance = 0;
end


%*************************************************************************
% 2. ���ٷ�֧������
%*************************************************************************
% ����֧����������Ч�ʡ�

% ע��: ���ҵĵ����ϣ�Matlab 2010b�汾�У�"for"����"vertcat"���Ч�ʸ��ߡ��Ҳ�֪��Ϊʲô��(LSSSSWC, 2011-07-25)
nViol   = zeros(N, 1);
violSum = zeros(N, 1);
for i = 1:N
    nViol(i)    = pop(i).nViol;
    violSum(i)  = pop(i).violSum;
end
% nViol   = vertcat(pop(:).nViol);
% violSum = vertcat(pop(:).violSum);

obj     = vertcat(pop(:).obj);
domMat  = calcDominationMatrix(nViol, violSum, obj); % �������Ч�ʵ�֧�����


% ����ÿ�������np��sp
for p = 1:N-1
    for q = p+1:N
        if(domMat(p, q) == 1)          % p ֧�� q
            ind(q).np = ind(q).np + 1;
            ind(p).sp = [ind(p).sp , q];
        elseif(domMat(p, q) == -1)     % q ֧�� p
            ind(p).np = ind(p).np + 1;
            ind(q).sp = [ind(q).sp , p];
        end
    end
end


% ��һ��ǰ�أ��ȼ� = 1��
front(1).f = [];    % 'front'�ṹ����ֻ��һ���ֶ�'f'��
                    % ���ǹ���ģ���Ϊÿ��ǰ���е�
                    % ���������ǲ�ͬ�ġ�
for i = 1:N
    if( ind(i).np == 0 )
        pop(i).rank = 1;
        front(1).f = [front(1).f, i];
    end
end

% ����ÿ������������еȼ����� pop(:).rank
fid = 1;        % ������ǰ��ID
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
front(fid) = [];    % ɾ�����һ���յ�ǰ�ؼ���



%*************************************************************************
% 3. �������
%*************************************************************************
if(isempty(opt.refPoints))
    pop = calcCrowdingDistance(opt, pop, front);
else
    [opt, pop] = calcPreferenceDistance(opt, pop, front);
end





function domMat = calcDominationMatrix(nViol, violSum, obj)
% ����: domMat = calcDominationMatrix(nViol, violSum, obj)
% ����: ����֧����󣬸þ���ʹ��Լ��֧���ϵָ����������
%   ֮���֧���ϵ��
%
% ����: 
%   domMat(N,N) : ֧�����
%       domMat(p,q)=1  : p ֧�� q
%       domMat(p,q)=-1 : q ֧�� p
%       domMat(p,q)=0  : ����֧��
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
        % 1. p �� q ���ǿ��е�
        %*************************************************************************
        if(nViol(p) == 0 && nViol(q)==0)
            pdomq = false;
            qdomp = false;
            for i = 1:numObj
                if( obj(p, i) < obj(q, i) )         % Ŀ�꺯������С����
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
        % 2. p ����, q ������
        %*************************************************************************
        elseif(nViol(p) == 0 && nViol(q)~=0)
            domMat(p, q) = 1;
        %*************************************************************************
        % 3. q ����, p ������
        %*************************************************************************
        elseif(nViol(p) ~= 0 && nViol(q)==0)
            domMat(p, q) = -1;
        %*************************************************************************
        % 4. p �� q ���ǲ����е�
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
% ����: [opt, pop] = calcPreferenceDistance(opt, pop, front)
% ����: ����R-NSGA-II��ʹ�õ�'ƫ�þ���'��
% ����: 
%   opt : ���� opt.refUseNormDistance=='ever' ʱ���˽ṹ����ܱ��޸ġ�
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.1  Data: 2011-07-26
%*************************************************************************

%*************************************************************************
% 1. ��ʼ��
%*************************************************************************
numObj = length( pop(1).obj );  % Ŀ������

refPoints = opt.refPoints;
refWeight = opt.refWeight;      % Ŀ��Ȩ������
if(isempty(refWeight))
    refWeight = ones(1, numObj);
end
epsilon = opt.refEpsilon;
numRefPoint = size(refPoints, 1);

% ������һ������
bUseFrontMaxMin = false;    % bUseFrontMaxMin : �Ƿ�ʹ��ǰ���е�������Сֵ��Ϊ��һ�����ӡ�
if( strcmpi(opt.refUseNormDistance, 'ever') )
    % 1) �ҵ�ÿ��Ŀ����ܵ�������Сֵ�������ڵ�ǰ��Ⱥ��
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
    % 2) ��ʹ�ù�һ����ŷ�Ͼ��롣
    bUseFrontMaxMin = true;
elseif( strcmpi(opt.refUseNormDistance, 'no') )
    % 3) ��ʹ�ù�һ����ŷ�Ͼ��롣
    objMaxMin = ones(1,numObj);
else
    % 3) ����
    error('NSGA2:ParamError', ...
        '��֧�ֵĲ���: options.refUseNormDistance="%s", ֻ֧�� "yes" �� "no"',...
        opt.refUseNormDistance);
end


%*************************************************************************
% 2. ����ƫ�þ��� pop(:).prefDistance
%*************************************************************************
for fid = 1:length(front)
    % ����1: ����ÿ��ǰ���еļ�Ȩŷ�Ͼ���
    idxFront = front(fid).f;            % idxFront : ��ǰǰ���и��������
    numInd = length(idxFront);          % numInd : ��ǰǰ���еĸ�������
    popFront = pop(idxFront);           % popFront : ��ǰ��fid�еĸ���

    objFront = vertcat(popFront.obj);   % objFront : ���и����ȫ��Ŀ��ֵ

    if(bUseFrontMaxMin)
        objMaxMin = max(objFront) - min(objFront); % objMaxMin : ��ǰǰ���еĹ�һ������
    end

    % normDistance : ��Ȩ��һ��ŷ�Ͼ���
    normDistance = calcWeightNormDistance(objFront, refPoints, objMaxMin, refWeight);
    
    
    % ����2: ����ƫ�þ���
    prefDistanceMat = zeros(numInd, numRefPoint);
    for ipt = 1:numRefPoint
        [~,ix] = sort(normDistance(:, ipt));
        prefDistanceMat(ix, ipt) = 1:numInd;
    end
    prefDistance = min(prefDistanceMat, [], 2);
    clear ix

    
    % ����3: Epsilon �������
    idxRemain = 1:numInd;           % idxRemain : δ������������
    while(~isempty(idxRemain))
        % 1. ��ʣ�������ѡ��һ��
        objRemain = objFront( idxRemain, :);
        selIdx = randi( [1,length(idxRemain)] );
        selObj = objRemain(selIdx, :);

        % 2. �����һ��ŷ�Ͼ���
        % distanceToSel : ����ѡ��Ĺ�һ��ŷ�Ͼ���
        distanceToSel = calcWeightNormDistance(objRemain, selObj, objMaxMin, refWeight);
        

        % 3. ����epsilon�����ڵĸ���
        idx = find( distanceToSel <= epsilon );     % idx : �� idxRemain �е�����
        if(length(idx) == 1)    % Ψһ�ĸ�����Ǳ�ѡ�е��Ǹ�
            idxRemain(selIdx)=[];
        else
            for i=1:length(idx)
                if( idx(i)~=selIdx )
                    idInIdxRemain = idx(i);     % idx ���� idxRemain �����е�����
                    id = idxRemain(idInIdxRemain);
                    
                    % *����ƫ�þ����Խ�����Щ����
                    % ��ѡ��ļ��ʡ�
                    prefDistance(id) = prefDistance(id) + round(numInd/2);
                end
            end
            idxRemain(idx) = [];
        end
        
    end

    % ����ƫ�þ���
    for i=1:numInd
        id = idxFront(i);
        pop(id).prefDistance = prefDistance(i);
    end
end


function distance = calcWeightNormDistance(points, refPoints, maxMin, weight)
% ����: calcWeightNormDistance(points, refPoints, maxMin, weight)
% ����: �����"points"��"refPoints"�ļ�Ȩŷ�Ͼ���
% ����: 
%   points(nPoint, N)       : ÿһ����Nά�ռ��е�һ���㡣
%   refPoints(nRefPoint, N) : ÿһ����һ���ο��㡣
%   maxMin(1, N)            : ��һ�����ӡ�
%   weight(1, N)            : Ȩ��
%
% ����: 
%   distance(nPoint, nRefPoint)
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-14
%*************************************************************************

nRefPoint = size(refPoints, 1);     % �ο�������
nPoint = size(points, 1);           % �������

distance = zeros(nPoint, nRefPoint);
for ipt = 1:nRefPoint
    refpt = refPoints(ipt, :);
    for i = 1:nPoint
        weightNormDist = ((points(i, :)-refpt) ./ maxMin).^2 .* weight;
        distance(i, ipt) = sqrt(sum(weightNormDist));
    end
end





function pop = calcCrowdingDistance(opt, pop, front)
% ����: pop = calcCrowdingDistance(opt, pop, front)
% ����: ����ԭʼNSGA-II��ʹ�õ�'ӵ���Ⱦ���'��
% �﷨:
% ����: 
% ����: 
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-11
%*************************************************************************

numObj = length( pop(1).obj );  % Ŀ������
for fid = 1:length(front)
    idx = front(fid).f;
    frontPop = pop(idx);        % frontPop : ��ǰ��fid�еĸ���
    
    numInd = length(idx);       % nInd : ��ǰǰ���еĸ�������
    
    obj = vertcat(frontPop.obj);
    obj = [obj, idx'];          % Ŀ��ֵ�����IDһ������
    for m = 1:numObj
        obj = sortrows(obj, m);

        colIdx = numObj+1;
        pop( obj(1, colIdx) ).distance = Inf;         % ��һ������
        pop( obj(numInd, colIdx) ).distance = Inf;    % ���һ������
        
        minobj = obj(1, m);         % Ŀ��m����Сֵ
        maxobj = obj(numInd, m);    % Ŀ��m�����ֵ
        
        for i = 2:(numInd-1)
            id = obj(i, colIdx);
            pop(id).distance = pop(id).distance + (obj(i+1, m) - obj(i-1, m)) / (maxobj - minobj);
        end
    end
end