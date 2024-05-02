function gwf_plot_analysis(gwf, rf, dt, opt)
% function gwf_plot_analysis(gwf, rf, dt, col)
%
% gwf - n x (1, 2, or 3)  [gradient as executed by amplifier]
% rf  - n x 1             [effect of RF pulses, either 1 or -1]
% dt  - 1 x 1             [time step]
% opt - options structure

if (nargin < 4), opt = []; end, opt = gwf_opt(opt);

txt = gwf_analysis(gwf, rf, dt, opt);

for c = 1:numel(txt)
    txt{c} = [txt{c} '^ _ '];
end

text(0,-5, txt); 
ylim([-20 10]);
axis off;








