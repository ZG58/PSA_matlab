function result = loadpopfile(fileName)
% 函数: result = loadpopfile(fileName)
% 描述: 加载上次优化生成的种群文件。
% 语法:
%       oldresult = loadpopfile('populations.txt');
% 返回: 
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-01
%*************************************************************************


%*************************************************************************
% 打开种群文件
%*************************************************************************
if( ~ischar(fileName))
    error('NSGA2:LoadPopError', '参数 fileName 应该为字符串！');
end

fid = fopen(fileName, 'r');
if(fid==-1)
    error('NSGA2:LoadPopError', '种群文件 "%s" 无法打开！', fileName);
end
%fprintf('...正在加载种群文件 "%s"......\n', fileName);



%*************************************************************************
% 读取文件头
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
        '种群文件 "%s" 不是一个Nsga2种群文件！行：\n%s\n', ...
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
            warning('NSGA2:PopFileError', '不支持的状态关键字: "%s"', token{1}{1});
    end
    
    strLine = fgetl(fid);
end


%*************************************************************************
% 初始化结果结构体
%*************************************************************************
pop = repmat( struct(...
    'var', zeros(1,nVar), ...
    'obj', zeros(1,nObj), ...
    'cons', zeros(1,nCons)),...
    [1,popsize]);

% state: 某一代替换状态
state = struct();
for i = 1:length(fieldNames)
    state.(fieldNames{i}) = 0;
end

result.pops     = repmat(pop, [maxGen, 1]);     % 每一行是某一代的种群
result.states   = repmat(state, [maxGen, 1]);   % 每一行是某一代的优化状态
clear i fieldNames pop


%*************************************************************************
% 解析种群文件
%*************************************************************************
lastGen = 0;    % 最后一个拥有正确数据的代数。
try
    strLine = fgetl(fid);
    while ischar(strLine)
        %*****************************************************************
        % 1. 跳过空行
        strLine = strtrim(strLine);
        if( isempty(strLine)  )
            strLine = fgetl(fid);   % 读取新行
            continue;
        end
        
        %*****************************************************************
        % 2. 某一代的第一行
        if( strcmp(strLine(1:11), '#Generation') == 1)
            % 只需要 'ngen'
            ngen = sscanf(strLine(12:end), ' %d');
        else
            error('NSGA2:PopFileError', '种群文件格式错误！行：\n%s\n', strLine);
        end


        %*****************************************************************
        % 3. 读取优化状态
        strLine = fgetl(fid);
        while( ischar(strLine) && strcmp(strLine, '#end')==0 )
            token = textscan(strLine, '%s%f');
            result.states(ngen).(token{1}{1}) = token{1,2};
            
            strLine = fgetl(fid);
        end
            
        
        %*****************************************************************
        % 4. 读取种群
        strLine = fgetl(fid);
        val = fscanf(fid, '%f');
        ncols = nVar+nObj+nCons;
        nrows = popsize;
        if( length(val) ~= ncols*nrows )
            error('NSGA2:PopFileError', '读取种群数据时文件出错！');
        end

        val = reshape(val, ncols, nrows)';  % 按列重塑向量
        for i = 1:popsize
            result.pops(ngen, i).var = val(i, 1:nVar);
            result.pops(ngen, i).obj = val(i, (nVar+1):(nVar+nObj));
            result.pops(ngen, i).cons= val(i, (nVar+nObj+1):end);
        end
        
        %*****************************************************************
        % 读取下一行
        lastGen = ngen;
        strLine = fgetl(fid);

    end
catch exception
    id = exception.identifier;
    msg = [exception.message, 'File = ', fileName];
    warning(id, msg);       % 如果文件末尾有误，返回正确的数据。
end


%*************************************************************************
% 进行一些清理工作
%*************************************************************************

% 删除未使用的数据
result.pops(lastGen+1:end, :) = [];
result.states(lastGen+1:end)  = [];
%fprintf('...加载种群文件成功。最后一代是 %d。\n', lastGen);


% 关闭文件
fclose(fid);