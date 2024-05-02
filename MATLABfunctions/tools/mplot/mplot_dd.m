function h = mplot_dd(D_list, w_list, color_list, dmin, dmax, h, opt)
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

if (nargin < 2), w_list = []; end
if (nargin < 3), color_list = []; end
if (nargin < 4), dmin = []; end
if (nargin < 5), dmax = []; end
if (nargin < 6), h = gca; end
if (nargin < 7), opt = []; end

% Autoconfigure options
opt = mplot_opt(opt);

a = numel(color_list)/3;

cla(h); hold(h, 'on');

% Final plot
for c = 1:a
    w = w_list{c} / sum(w_list{c});
    D = log10(D_list{c});
    [D_sorted, D_order] = sort(D);
    w_sorted = w(D_order);
    
    stem(h,D_sorted,w_sorted,'-o', 'filled', 'Color', color_list(c,:), 'LineWidth',1, 'MarkerSize',5);
    
end
    
set(h,...
    'XLim',[log10(dmin)-.1 log10(dmax)+.1], ...
    'YLim',[0 1], ...
    'YAxisLocation','left', ...
    'XTick', -12:.5:-7, ...
    'YTick',0:.1:1, ...
    'TickDir','out', ...
    'TickLength',.02*[1 1], ...
    'FontSize',opt.mplot.fs, ...
    'LineWidth',opt.mplot.lw, ...
    'Box','on')

axis(h,'square');

xlabel(h,'\rmlog_1_0(\itD\rm / m^2s^-^1)','FontSize',opt.mplot.fs)
ylabel(h,'\it{P}(\it{D})','FontSize',opt.mplot.fs)

