function dtd_covariance_plot(S, xps, h, h2)
% function dtd_covariance_plot(S, xps, h, h2)

if (nargin < 3), h  = gca; end
if (nargin < 4), h2 = []; end

opt = dtd_covariance_opt();
opt = mplot_opt(opt);
opt = mdm_opt(opt);

[m,~,n] = dtd_covariance_1d_data2fit(S, xps, opt);
S_fit = dtd_covariance_1d_fit2data(m, xps);

% s0
s0 = m(1);

% diffusion tensor
dt_1x6 = m(2:7)' * 1e9;
dps = tm_dt_to_dps(dt_1x6);

% covariance tensor
ct_1x21  = m(8:28)' * 1e18;
dps = tm_ct_to_dps(ct_1x21, dps);

% plot with powder average, but also visualize the directions
xps_pa = mdm_xps_pa(msf_rmfield(xps, 'c_volume'));
S_pa = squeeze(mio_pa(permute(cat(2, S, S_fit), [2 3 4 1]), xps));


c_case = 1;
switch (c_case)
    
    case 1       
        
        % find unique b_deltas
        b_delta_input = round(xps_pa.b_delta * 100) / 100;
        [b_delta_uni,~] = unique(b_delta_input);
        b_delta_uni = sort(b_delta_uni);
        
        
        % plotting color scheme
        if (numel(b_delta_uni) == 1)
            cmap = [0 0 0];
        else
            cmap = 0.7 * hsv(numel(b_delta_uni) + 2);
        end
                       
                
        for c = 1:numel(b_delta_uni)
            
            % get all data for this b_delta
            ind = abs(xps_pa.b_delta - b_delta_uni(c)) < opt.mdm.pa.db_delta2;
            
            b = xps_pa.b(ind);
            
            semilogy(h, b * 1e-9, S_pa(1,ind) / s0, 'o', ...
                'markersize', 5, 'color', cmap(c,:), ...
                'markerfacecolor', cmap(c,:)); 
            hold(h, 'on');
            
            % get a nice and smooth prediction
            p = polyfit(xps_pa.b(ind), log(S_pa(2,ind))', 2);
            b2 = linspace(0, max(xps.b(:)) * 1.05, 100);
            s2 = exp(polyval(p, b2));
            
            semilogy(h, b2 * 1e-9, s2 / s0, '-', 'linewidth', 2, 'color', cmap(c,:));
  
            
            % plot the actual datapoints
            for c_b = 1:numel(b)
                
                ind2 = abs(xps.b_delta - b_delta_uni(c)) < opt.mdm.pa.db_delta2;
                ind2 = ind2 & abs(xps.b - b(c_b)) < opt.mdm.pa.db;
                

                % see how well aligned it is with the primary direction
                p = linspace(-1, 1, sum(ind2))';
                p = p / 10;
                
                [~,ind3] = sort(S_fit(ind2));
                sf = @(x) x(ind3);
                
                
                plot(h, (b(c_b) * 1e-9 + p), sf(S(ind2)) / s0, 'k.', 'color', cmap(c,:));
                plot(h, (b(c_b) * 1e-9 + p), sf(S_fit(ind2)) / s0, 'k-', 'color', cmap(c,:));
                
            end
            
            
        end
        
        xlim(h, [0 max(xps.b) * 1.1e-9]);
        ylim(h, [0.05 1.1]);        
    
        
        set(h,...
            'TickDir','out', ...
            'TickLength',.02*[1 1], ...
            'FontSize',opt.mplot.fs, ...
            'LineWidth',opt.mplot.lw, ...
            'Box','off')
        
    
    case 2 % old plot
        
        plot(h, S, 'k.');
        hold(h, 'on');
        plot(h, S_fit, 'r-');
        xlim(h, [0 1+numel(S_fit)]);
        ylim(h, [0 max(max(S_fit), max(S)) * 1.1]);
end


title(h, sprintf('MD = %1.2f, FA = %1.2f\nMK_I = %1.2f, MK_A = %1.2f\n  MEAS RANK = %i', ...
    dps.MD, dps.FA, dps.MKi, dps.MKa, n));


if ~isempty(h2)
    
    C = m(8:end) * 1e18;
    
    [x,y,z] = sphere(120);
    
    c = zeros(size(x));
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            r = [x(i,j) y(i,j) z(i,j)];
            
            r4 = tm_1x6_to_1x21(tm_1x3_to_1x6(1,0,r));
            r4 = r4 / sqrt(tm_inner(r4, r4));
            L1 = tm_inner(C, r4);
            
            L = L1;
            
            x(i,j) = x(i,j) * L;
            y(i,j) = y(i,j) * L;
            z(i,j) = z(i,j) * L;
            
            c(i,j) = L;
            
        end
    end
    
    surface(h2,x,z,y,c);
    
    caxis(h2, [0 1.5] );
    colormap(h2, fliplr(jet(100)));
    axis(h2, 'off')
    axis(h2, 'tight');
    axis(h2, 'equal');
    axis(h2, 'vis3d');
    camlight(h2, 0, 0);

    material([0.7 0.4 0.8 1]); 
    shading(h2, 'interp');

end