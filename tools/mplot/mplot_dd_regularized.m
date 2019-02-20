function h = mplot_dd_regularized(D, PD_list, color_list, h, opt)
% function h = mplot_dtd(t_1x6, w, dmin, dmax, h)
%
% Plot the distribution of diffusion tensors in a 2D plot
%
% t_1x6    - list of discrete tensors
% w        - relative weight of each tensor in the distribution
% dmin     - min isotropic diffusivity (for determining axis limits)
% dmax     - max isotropic diffusivity
% h        - reference to the plot area
% opt      - options
%            for powder averaged plotting, set opt.mplot.dtd_col_mode = 0

if (nargin < 2), PD_list = []; end
if (nargin < 3), color_list = []; end
if (nargin < 4), h = gca; end
if (nargin < 5), opt = []; end

% Autoconfigure options
opt = mplot_opt(opt);

a = numel(color_list)/3;

cla(h); hold(h, 'on');

log_D = log10(D);
[D_sorted, D_order] = sort(log_D);
log_dmin = round(min(log_D),1)-0.1;
log_dmax = round(max(log_D),1)+0.1;

max_PD = 0;

% Final plot
for c = 1:a
    PD_sorted = PD_list{c}(D_order);
    if max(PD_sorted) > max_PD
        max_PD = max(PD_sorted);
    end
    plot(h,D_sorted,PD_sorted,'-', 'Color', color_list(c,:), 'LineWidth',2);
end
    
max_PD = round(max_PD+0.01,2);

set(h,...
    'XLim',[log_dmin log_dmax], ...
    'YLim',[0 max_PD], ...
    'YAxisLocation','left', ...
    'XTick', log_dmin:.5:log_dmax, ...
    'YTick',0:.05:max_PD, ...
    'TickDir','out', ...
    'TickLength',.02*[1 1], ...
    'FontSize',opt.mplot.fs, ...
    'LineWidth',opt.mplot.lw, ...
    'Box','on')

axis(h,'square');

xlabel(h,'\rmlog_1_0(\itD\rm / m^2s^-^1)','FontSize',opt.mplot.fs)
ylabel(h,'\it{P}(\it{D})','FontSize',opt.mplot.fs)

