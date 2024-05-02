function opt = mplot_opt(opt)
% function opt = mplot_opt(opt)
%
% Specifies default options 

if ((nargin < 1) || (isempty(opt)) || (~isfield(opt,'mplot')))
    opt.mplot.present = 1; 
end

opt.mplot = msf_ensure_field(opt.mplot, 'lw', 1);
opt.mplot = msf_ensure_field(opt.mplot, 'fs', 10);
opt.mplot = msf_ensure_field(opt.mplot, 'ms', 5);

opt.mplot = msf_ensure_field(opt.mplot, 'dtd_col_mode', 1);
opt.mplot = msf_ensure_field(opt.mplot, 'dtd_plot_type', 'point_estimate');

opt.mplot = msf_ensure_field(opt.mplot, 'terminology', 'topgaard17'); %lasic14, szczepankiewicz16, westin16, topgaard17




