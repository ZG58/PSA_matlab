function defaultopt = nsgaopt()
% ����: defaultopt = nsgaopt()
% ����: ����NSGA-IIĬ��ѡ��ṹ�塣
% �﷨:  opt = nsgaopt()
%         LSSSSWC, NWPU
%   Revision: 1.3  Data: 2011-07-13
%*************************************************************************


defaultopt = struct(...
... % �Ż�ģ��
    'popsize', 50,...           % ��Ⱥ��С
    'maxGen', 100,...           % ������
    'numVar', 0,...             % ��Ʊ�������
    'numObj', 0,...             % Ŀ������
    'numCons', 0,...            % Լ������
    'lb', [],...                % ��Ʊ����½� [1:numVar]
    'ub', [],...                % ��Ʊ����Ͻ� [1:numVar]
    'vartype', [],...           % ������������ [1:numVar]��1=ʵ��, 2=����
    'objfun', @objfun,...       % Ŀ�꺯��
... % �Ż�ģ�͸��������
    'nameObj',{{}},...
    'nameVar',{{}},...
    'nameCons',{{}},...
... % ��ʼ�������
    'initfun', {{@initpop}},...         % ��Ⱥ��ʼ������ (Ĭ��Ϊ�����)
    'outputfuns',{{@output2file}},...   % �������
    'outputfile', 'populations.txt',... % ����ļ���
    'outputInterval', 1,...             % ������
    'plotInterval', 5,...               % ���ε��� "plotnsga" ֮��ļ����
... % �Ŵ��㷨����
    'crossover', {{'intermediate', 1.2}},...         % �������� (����=1.2)
    'mutation', {{'gaussian',0.1, 0.5}},...          % �������� (��������=0.1, ��������=0.5)
    'crossoverFraction', 'auto', ...                 % �����б����Ľ�����
    'mutationFraction', 'auto',...                   % �����б����ı�����
... % �㷨����
    'useParallel', 'no',...                          % ���м�����Ⱥ��Ŀ�꺯���� {'yes','no'}
    'poolsize', 0,...                                % ���м���ʹ�õĹ�����������0 = �Զ�ѡ��
... % R-NSGA-II ����
    'refPoints', [],...                              % ����ָ��ƫ�õĲο��㡣ÿһ����һ���ο��㡣
    'refWeight', [],...                              % ����ŷ�Ͼ���ʱʹ�õ�Ȩ������
    'refUseNormDistance', 'front',...                % ʹ���ɿ��ܵ�������СĿ��ֵ��һ����ŷ�Ͼ��롣 {'front','ever','no'}
    'refEpsilon', 0.001 ...                          % ����epsilon��ѡ�������ʹ�õĲ���
);