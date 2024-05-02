function msf_fprintf(opt, str, varargin)
% function msf_fprintf(opt, str, varargin)
%
% prints to terminal if opt.verbose = 1

if (opt.verbose == 1)
    fprintf(str, varargin{:});
end
