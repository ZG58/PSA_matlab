function opt = output2file(opt, state, pop, type, varargin)
% 函数: opt = output2file(opt, state, pop, type, varargin)
% 描述: 将种群'pop'输出到文件。文件名由
% 'opt.outputfile' 字段指定。
% 参数: 
%   type : 输出类型。 -1 = 最后一次调用，关闭已打开的文件。
%          其他值(或不存在) = 正常输出
%   varargin : 在 options.outputfuns 元胞数组中定义的任何参数。
%
%         LSSSSWC, NWPU
%    Revision: 1.2  Data: 2011-07-13
%*************************************************************************


if(isempty(opt.outputfile))
    return;  % 未指定输出文件名，直接返回
end

if( isfield(opt, 'outputfileFID') )
    fid = opt.outputfileFID;
else
    fid = [];
end

%*************************************************************************
% 1.打开输出文件并输出一些种群信息
%*************************************************************************
if( isempty(fid) )
    fid = fopen(opt.outputfile, 'w');
    if( fid == 0)
        error('NSGA2:OutputFileError', '无法打开输出文件!! 文件名:%s', opt.outputfile);
    end
    opt.outputfileFID = fid;
    
    % 输出一些信息
    fprintf(fid, '#NSGA2\r\n');

    fprintf(fid, 'popsize %d\r\n', opt.popsize);
    fprintf(fid, 'maxGen %d\r\n', opt.maxGen);
    fprintf(fid, 'numVar %d\r\n', opt.numVar);
    fprintf(fid, 'numObj %d\r\n', opt.numObj);
    fprintf(fid, 'numCons %d\r\n', opt.numCons);
    
    % 输出状态字段名称
    fprintf(fid, 'stateFieldNames\t');
    names = fieldnames(state);
    for i = 1:length(names)
        fprintf(fid, '%s\t', names{i});
    end
    fprintf(fid, '\r\n');
    
    fprintf(fid, '#end\r\n\r\n\r\n');
end

%*************************************************************************
% 2. 如果这是最后一次调用，关闭输出文件
%*************************************************************************
if(type == -1)
    fclose(fid);
    rmfield(opt, 'outputfileFID');
    return
end

%*************************************************************************
% 3. 将种群输出到文件
%*************************************************************************
fprintf(fid, '#Generation %d / %d\r\n', state.currentGen, opt.maxGen);

% 输出每个状态字段
names = fieldnames(state);
for i = 1:length(names)
    fprintf(fid, '%s\t%g\r\n', names{i}, getfield(state, names{i}));
end
fprintf(fid, '#end\r\n');

for i = 1:opt.numVar
    fprintf(fid, 'Var%d\t', i);
end
for i = 1:opt.numObj
    fprintf(fid, 'Obj%d\t', i);
end
for i = 1:opt.numCons
    fprintf(fid, 'Cons%d\t', i);
end
fprintf(fid, '\r\n');

for p = 1 : opt.popsize
    for i = 1:opt.numVar
        fprintf(fid, '%g\t', pop(p).var(i) );
    end
    for i = 1:opt.numObj
        fprintf(fid, '%g\t', pop(p).obj(i) );
    end
    for i = 1:opt.numCons
        fprintf(fid, '%g\t', pop(p).cons(i));
    end
    fprintf(fid, '\r\n');
end

fprintf(fid, '\r\n\r\n\r\n');