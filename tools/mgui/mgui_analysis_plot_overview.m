function mgui_analysis_plot_overview(S, xps, h_top, h_bottom, S_fit)
% function mgui_analysis_plot_overview(S, xps, h_top, h_bottom)

if (isempty(S))
    axes(h_top);
    cla;
    text(0, 0, 'Draw ROI');
    ylim([-1 1]);
    xlim([0 10]);
    axis(h_top, 'off');
    axis(h_bottom, 'off');
    return;
end

col_gray = [0.0 0.0 0.0] + 0.7;
col_red  = [0.8 0.0 0.0];

xps = msf_ensure_field(xps, 'n', size(S,2));

x = 1:xps.n;
MS = mean(S, 1)';

hold(h_top, 'on');

if (size(S, 1) < 50)
    
    % Deal with complex signals
    if (any(~isreal(S(:))))
        S_plot = abs(S');
    else
        S_plot = S';
    end
    
    plot(h_top, x, S_plot, '.-', 'color', col_gray);
    
else
    
    if (any(~isreal(S(:))))
        MS_plot = abs(MS);
        MD_plot = std(abs(S), [], 1);
    else
        MS_plot = MS;
        MD_plot = std(S, [], 1);
    end
    
    errorbar(h_top, x, MS_plot, MD_plot, 'ks-');
    
    
    if (isfield(xps, 'c_volume'))
        plot(h_top, x(xps.c_volume), MS_plot(xps.c_volume), 'ko', 'markerfacecolor', 'red', 'markersize', 9);
    end
end

if (nargin > 4)
    plot(h_top, x, S_fit, '-', 'color', col_red);
end


axis(h_top, 'on');
box(h_top, 'off');
set(h_top, 'tickdir','out', 'ticklength', [0.03 0.1]);

% Compute y-scaling
if (any(~isreal(S(:))))
    y_axis = [0 max(abs(S(:)) + eps) * 1.1];
elseif (min(S(:)) < 0)
    y_axis = [-1 1] * (max(abs(S(:))) + eps) * 1.1;
else % standard case
    y_axis = [0 max(S(:) + eps) * 1.1];
end

axis(h_top, [ [1 xps.n] + [-1 1] * 0.5 y_axis]);

xlabel(h_top, 'Acq number');
ylabel(h_top, 'Signal');
set(h_top, 'xtick', [1 xps.n]);

f = @(x) num2str(getfield(msf_ensure_field(xps, x, NaN), x));

str = {'Info:', ...
    sprintf('xps.n: %s', f('n'))};

text(0,0,str, 'parent', h_bottom);
axis(h_bottom, 'off');


