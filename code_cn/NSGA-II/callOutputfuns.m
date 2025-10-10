function opt = callOutputfuns(opt, state, pop, type)
% 函数: opt = callOutputfuns(opt, state, pop, type)
% 描述: 调用输出函数(如果存在)。
% 参数: 
%   type : 输出类型。
%       -1 = 最后一次调用 (例如关闭文件)
%       其他值(或不存在) = 正常输出
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-13
%*************************************************************************


if(nargin <= 3)
    type = 0;   % 正常输出
end


if( ~isempty(opt.outputfuns) )
    fun = opt.outputfuns{1};
    opt = fun(opt, state, pop, type, opt.outputfuns{2:end});
end