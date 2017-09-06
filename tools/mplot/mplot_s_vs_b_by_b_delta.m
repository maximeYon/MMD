function m = mplot_s_vs_b_by_b_delta(S, xps, fun_data2fit, fun_fit2data, fun_opt, h, h2, ind)
% function mplot_s_vs_b_by_b_delta(S, xps, fun_data2fit, fun_fit2data, fun_opt, h, h2)

if (nargin < 3), fun_data2fit = []; end
if (nargin < 4), fun_fit2data = []; end
if (nargin < 5), fun_opt = []; end
if (nargin < 6), h = gca; end
if (nargin < 7), h2 = []; end
if (nargin < 8), ind = []; end

% adapt xps
xps = msf_ensure_field(xps, 'b_eta', zeros(size(xps.b_delta)));

% fit function

if (nargin < 8)   
    ind = ones(size(S)) == 1;
end

if (~isempty(fun_data2fit))
    m = fun_data2fit(S(ind), mdm_xps_subsample(xps, ind), fun_opt());
else
    m = [];
end



% find unique b_delta
b_delta_input = round(xps.b_delta * 100) / 100;

[b_delta_uni,~] = unique(b_delta_input);

b_delta_uni = sort(b_delta_uni);

if (numel(b_delta_uni) == 1)
    cmap = [0 0 0];
else
    cmap = 0.8 * hsv(numel(b_delta_uni) + 2);
end

cla(h);
hold(h, 'off');
for c = 1:numel(b_delta_uni)
    
    ind2 = b_delta_input == b_delta_uni(c);
    
    ind3 = ind2 & (~ind);
    
    % Construct a fitted signal
    clear xps2;
    try
        xps2.n = 32;
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
    end
    
    if (isempty(m))
        S_fit = zeros(xps2.n,1);
        sc = max(S(:));
    else
        try
            
            S_fit = fun_fit2data(m, xps2);
            
            % This now allows for separate scaling of each series
            sc = S_fit(1);
        catch
            xps2 = mdm_xps_subsample(xps, ind2);
            S_fit = fun_fit2data(m, xps2);
            
            % Scale by m(1) assuming it represents S0
            sc = m(1);
        end
    end
    
    semilogy(h, xps.b(ind2) * 1e-9, S(ind2) / sc, 'o', 'color', cmap(c,:), ...
        'markersize', 5, 'markerfacecolor', cmap(c,:));

    hold(h, 'on');

    semilogy(h, xps.b(ind3) * 1e-9, S(ind3) / sc, 'x', 'color', cmap(c,:), ...
        'markersize', 10, 'markerfacecolor', cmap(c,:));

    semilogy(h, xps2.b * 1e-9, S_fit / sc, '-', 'color', cmap(c,:), 'linewidth', 2);
end

axis(h, 'on');
box(h, 'off');
set(h, 'tickdir','out');


xlim(h, [0 max(xps2.b) * 1.1e-9]);
ylim(h, [10^floor(log10(min(S(S(:) > 0) / sc))) 1.1]);

xlabel(h, 'b [ms/um^2]');
ylabel(h, 'Signal');



% show legend
if (~isempty(h2))
    
    l_str = cell(1,numel(b_delta_uni));
    for c = 1:numel(b_delta_uni)
        plot(h2,10,10,'-', 'color', cmap(c,:), 'linewidth', 2);
        hold(h2, 'on');
        l_str{c} = ['b_\Delta = ' num2str(b_delta_uni(c))];
    end
    ylim(h2, [-1 1]);
    legend(h2, l_str, 'location', 'northwest');
    axis(h2, 'off');

end