function h = mplot_dtd(t_1x6, w, dmin, dmax, h, opt)
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

if (nargin < 2), w = []; end
if (nargin < 3), dmin = []; end
if (nargin < 4), dmax = []; end
if (nargin < 5), h = gca; end
if (nargin < 6), opt = []; end

% Autoconfigure weights: assume equal weights
if (isempty(w)), w = ones(size(t_1x6,1), 1); end

% Autoconfigure axis limits if needed
if (isempty(dmin) || isempty(dmax))
    iso = t_1x6 * [1 1 1 0 0 0] / 3;
    dmin = round(min(iso) * 1e9 * 10) / 10;
    dmax = round(max(iso) * 1e9 * 10) / 10;
end

% Autoconfigure options
opt = mplot_opt(opt);



ratiomax = dmax/dmin;

xmin = log10(dmin);
xmax = log10(dmax);
ymin = log10(1/ratiomax);
ymax = log10(ratiomax);

% convert weights to marker sizes
w = w / sum(w);
ms = 5 * opt.mplot.ms * sqrt(w);

if (max(ms) > 15), ms = ms / 15 * 15; end

cla(h); hold(h, 'on');
for c = 1:size(t_1x6) % loop over the tensors in the distribution
    
    p = tm_3x3_to_tpars(tm_1x6_to_3x3(t_1x6(c,:)));
    
    x = log10(p.iso);
    y = log10(p.para/p.perp);
    
    if (opt.mplot.dtd_col_mode > 0)
        col = p.fa * abs(p.norm);
    else
        col = 'black';
    end
    
    plot(h,x,y,'o',...
        'LineWidth',.01, ...
        'MarkerSize',ms(c),...
        'Color',col,...
        'MarkerFaceColor',col);
end

set(h,...
    'XLim',[xmin xmax], ...
    'YLim',[ymin ymax], ...
    'YAxisLocation','right', ...
    'XTick', -11:.5:-8, ...
    'YTick',-2:.5:2, ...
    'TickDir','out', ...
    'TickLength',.02*[1 1], ...
    'FontSize',opt.mplot.fs, ...
    'LineWidth',opt.mplot.lw, ...
    'Box','on')

axis(h,'square');

xlabel(h,'size, log \it{D}\rm_{iso} / m^2s^-^1','FontSize',opt.mplot.fs)
ylabel(h,'shape, log \it{D}\rm_{para} / \it{D}\rm_{perp}','FontSize',opt.mplot.fs)

if (opt.mplot.dtd_col_mode > 0)
    title(h,'orientation, [RGB]=[xyz]_ ','FontSize',opt.mplot.fs,'FontWeight','normal')
end

