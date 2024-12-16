function mplot_s_vs_b_by_b_delta_ilt_regularized(S, xps, fun_data2fit, fun_fit2data, opt, h, h2, ind)
% function mplot_s_vs_b_by_b_delta(S, xps, fun_data2fit, fun_fit2data, opt, h, h2)

if (nargin < 3), fun_data2fit = []; end
if (nargin < 4), fun_fit2data = []; end
if (nargin < 5), opt = []; end
if (nargin < 6), h = gca; end
if (nargin < 7), h2 = []; end
if (nargin < 8), ind = []; end

% Hack to make things work again
if (isa(opt, 'function_handle'))
    opt = opt();
end

% init
if (isa(opt, 'function_handle')), opt = opt(); end
opt = mplot_opt(opt);

data_marker = 'o';

if (isempty(fun_data2fit) || isempty(fun_fit2data))
    data_marker = 'o-';
end

% adapt xps
xps = msf_ensure_field(xps, 'b_eta', zeros(size(xps.b_delta)));

% This is for powder averaged data... so make it powder averaged if needed
xps_pa = mdm_xps_pa(msf_rmfield(xps, 'c_volume'));

if (xps_pa.n ~= xps.n)
    if (isfield(xps,'c_volume'))
        c_volume = xps.c_volume;
        xps = rmfield(xps, 'c_volume');
    else
        c_volume = [];
    end
    
    [S,xps] = mdm_powder_average_1d(S, xps);
    
    xps.c_volume = c_volume;
end

% fit function
if (nargin < 8)
    ind = ones(size(S)) == 1;
end

% find unique b_delta
b_delta_input = round(xps.b_delta * 100) / 100;
[b_delta_uni,~] = unique(b_delta_input);
b_delta_uni = sort(b_delta_uni);

if (numel(b_delta_uni) == 1)
    cmap = [0 0 0];
else
    cmap = 0.7 * hsv(numel(b_delta_uni) + 2);
end

% Allow 'hold on' to be used
try
    if (~strcmp(get(h, 'nextplot'), 'add'))
        error('jump out');
    end
catch
    cla(h);
    hold(h, 'off');
end

% Plot for legend
hold(h, 'off');
Nlines = numel(b_delta_uni);
for c = 1:Nlines
    lw = 2*(Nlines-.5*(c-1))/Nlines;
    semilogy(h, -1,2,'-', 'color', cmap(c,:), 'linewidth', lw); hold(h,'on');
end

m_list = { };

for c = 1:numel(b_delta_uni)
    
    ind2 = b_delta_input == b_delta_uni(c);
    
    if (~isempty(fun_data2fit))
        [D,dD,PD] = fun_data2fit(S(ind2), mdm_xps_subsample(xps, ind2), opt);
        m_list{c} = PD/sum(PD);
    else
        m = [];
        D = [];
        dD = [];
        m_list{c} = [];
        sqrt_chisq = [];
    end
    
    % Construct a fitted signal
    clear xps2;
    try
        xps2.n = 64;
        xps2.b = linspace(eps, max(xps.b(ind2)) * 1.1, xps2.n)';
        
        z = zeros(size(xps2.b));
        
        xps2.b_delta = z + mean(xps.b_delta(ind2));
        xps2.b_eta   = z + mean(xps.b_eta(ind2));
        
        if (isfield(xps, 's_ind'))
            tmp = xps.s_ind(ind2);
            
            if (numel(unique(tmp)) ~= 1)
                error('Cannot figure out which s_ind to use');
            end
            
            xps2.s_ind = z + tmp(1); % all should be equal
        end
    catch me
        1;
    end
    
    if (isempty(m_list{c}))
        S_fit = zeros(xps2.n,1)/max(S(:));
    else
        try
            S_fit = fun_fit2data(PD, xps2);
        catch
            xps2 = mdm_xps_subsample(xps, ind2);
            S_fit = fun_fit2data(PD, xps2);
        end
    end
    
    ms = 10*(Nlines-.5*(c-1))/Nlines;
    lw = 4*(Nlines-.9*(c-1))/Nlines;
    
    semilogy(h, xps.b(ind2) * 1e-9, S(ind2) / max(S(ind2)), data_marker, 'color', cmap(c,:), ...
        'markersize', ms, 'markerfacecolor', 'none', 'linewidth', 1);
    
    hold(h, 'on');
    
%     semilogy(h, xps.b(ind3) * 1e-9, S(ind3) / sc, 'x', 'color', cmap(c,:), ...
%         'markersize', 10, 'markerfacecolor', cmap(c,:));
    semilogy(h, xps2.b * 1e-9, S_fit, '-', 'color', cmap(c,:), 'linewidth', lw);
end

axis(h, 'on');

xlim(h, [0 max(xps2.b) * 1.1e-9]);
ylim(h, [.02 1.1]);
xlabel(h, ['\itb\rm [ms/',char(181),'m^2]'], 'fontsize', opt.mplot.fs);
ylabel(h, 'Signal', 'fontsize', opt.mplot.fs);

% show legend
l_str = cell(1,numel(b_delta_uni));
for c = 1:numel(b_delta_uni)
    plot(h,10,10,'-', 'color', cmap(c,:), 'linewidth', 2);
    hold(h, 'on');
    l_str{c} = ['\itb\rm_\Delta = ' num2str(b_delta_uni(c))];
end
legend(h, l_str, 'location', 'northeast');
legend(h, 'boxoff');

set(h,...
    'TickDir','out', ...
    'TickLength',.02*[1 1], ...
    'FontSize',opt.mplot.fs, ...
    'LineWidth',opt.mplot.lw, ...
    'Box','off')

name = char(fun_data2fit);
method_name = name(1:(strfind(name, '1d_data2fit')-2));

color_list = cmap(1:numel(b_delta_uni),:);

if (~isempty(m_list)) && strcmp(method_name,'ilt_regularized')
    mplot_dd_regularized(D, m_list, color_list, h2, opt);
    
else
    axis(h2, 'off');
end

end

