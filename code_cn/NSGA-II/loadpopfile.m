function result = loadpopfile(fileName)
% ����: result = loadpopfile(fileName)
% ����: �����ϴ��Ż����ɵ���Ⱥ�ļ���
% �﷨:
%       oldresult = loadpopfile('populations.txt');
% ����: 
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-01
%*************************************************************************


%*************************************************************************
% ����Ⱥ�ļ�
%*************************************************************************
if( ~ischar(fileName))
    error('NSGA2:LoadPopError', '���� fileName Ӧ��Ϊ�ַ�����');
end

fid = fopen(fileName, 'r');
if(fid==-1)
    error('NSGA2:LoadPopError', '��Ⱥ�ļ� "%s" �޷��򿪣�', fileName);
end
%fprintf('...���ڼ�����Ⱥ�ļ� "%s"......\n', fileName);



%*************************************************************************
% ��ȡ�ļ�ͷ
%*************************************************************************
popsize = 1;
maxGen  = 1;
nVar  = 0;
nObj  = 0;
nCons = 0;
fieldNames = {''};


strLine = fgetl(fid);
if( ~ischar(strLine) || strcmp(strLine, '#NSGA2')==0 )
    error('NSGA2:PopFileError', ...
        '��Ⱥ�ļ� "%s" ����һ��Nsga2��Ⱥ�ļ����У�\n%s\n', ...
        fileName, strLine);
end

strLine = fgetl(fid);
while( ischar(strLine) && strcmp(strLine, '#end')==0)
    token = textscan(strLine, '%s');
    keyword = strtrim(token{1}{1});
    switch keyword
        case 'popsize'
            popsize = str2double(token{1}{2});
        case 'maxGen'
            maxGen = str2double(token{1}{2});
        case 'numVar'
            nVar = str2double(token{1}{2});
        case 'numObj'
            nObj = str2double(token{1}{2});
        case 'numCons'
            nCons = str2double(token{1}{2});
        case 'stateFieldNames'
            nfield = length(token{1});
            fieldNames = cell(nfield - 1, 1);
            [fieldNames{:}] = token{1}{2:end};
        otherwise
            warning('NSGA2:PopFileError', '��֧�ֵ�״̬�ؼ���: "%s"', token{1}{1});
    end
    
    strLine = fgetl(fid);
end


%*************************************************************************
% ��ʼ������ṹ��
%*************************************************************************
pop = repmat( struct(...
    'var', zeros(1,nVar), ...
    'obj', zeros(1,nObj), ...
    'cons', zeros(1,nCons)),...
    [1,popsize]);

% state: ĳһ���滻״̬
state = struct();
for i = 1:length(fieldNames)
    state.(fieldNames{i}) = 0;
end

result.pops     = repmat(pop, [maxGen, 1]);     % ÿһ����ĳһ������Ⱥ
result.states   = repmat(state, [maxGen, 1]);   % ÿһ����ĳһ�����Ż�״̬
clear i fieldNames pop


%*************************************************************************
% ������Ⱥ�ļ�
%*************************************************************************
lastGen = 0;    % ���һ��ӵ����ȷ���ݵĴ�����
try
    strLine = fgetl(fid);
    while ischar(strLine)
        %*****************************************************************
        % 1. ��������
        strLine = strtrim(strLine);
        if( isempty(strLine)  )
            strLine = fgetl(fid);   % ��ȡ����
            continue;
        end
        
        %*****************************************************************
        % 2. ĳһ���ĵ�һ��
        if( strcmp(strLine(1:11), '#Generation') == 1)
            % ֻ��Ҫ 'ngen'
            ngen = sscanf(strLine(12:end), ' %d');
        else
            error('NSGA2:PopFileError', '��Ⱥ�ļ���ʽ�����У�\n%s\n', strLine);
        end


        %*****************************************************************
        % 3. ��ȡ�Ż�״̬
        strLine = fgetl(fid);
        while( ischar(strLine) && strcmp(strLine, '#end')==0 )
            token = textscan(strLine, '%s%f');
            result.states(ngen).(token{1}{1}) = token{1,2};
            
            strLine = fgetl(fid);
        end
            
        
        %*****************************************************************
        % 4. ��ȡ��Ⱥ
        strLine = fgetl(fid);
        val = fscanf(fid, '%f');
        ncols = nVar+nObj+nCons;
        nrows = popsize;
        if( length(val) ~= ncols*nrows )
            error('NSGA2:PopFileError', '��ȡ��Ⱥ����ʱ�ļ�����');
        end

        val = reshape(val, ncols, nrows)';  % ������������
        for i = 1:popsize
            result.pops(ngen, i).var = val(i, 1:nVar);
            result.pops(ngen, i).obj = val(i, (nVar+1):(nVar+nObj));
            result.pops(ngen, i).cons= val(i, (nVar+nObj+1):end);
        end
        
        %*****************************************************************
        % ��ȡ��һ��
        lastGen = ngen;
        strLine = fgetl(fid);

    end
catch exception
    id = exception.identifier;
    msg = [exception.message, 'File = ', fileName];
    warning(id, msg);       % ����ļ�ĩβ���󣬷�����ȷ�����ݡ�
end


%*************************************************************************
% ����һЩ������
%*************************************************************************

% ɾ��δʹ�õ�����
result.pops(lastGen+1:end, :) = [];
result.states(lastGen+1:end)  = [];
%fprintf('...������Ⱥ�ļ��ɹ������һ���� %d��\n', lastGen);


% �ر��ļ�
fclose(fid);