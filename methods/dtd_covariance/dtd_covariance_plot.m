function dtd_covariance_plot(S, xps, axh, axh2)
% function dtd_covariance_plot(S, xps, h, h2)

if (nargin < 4), axh2 = []; end

opt = dtd_covariance_opt();
%[m,~,n] = dtd_covariance_1d_data2fit(S, xps, opt);
%S_fit = dtd_covariance_1d_fit2data(m, xps);
m = mplot_signal_and_fit(S, xps, @dtd_covariance_1d_data2fit, @dtd_covariance_1d_fit2data, axh, opt);


% diffusion tensor
dt_1x6 = m(2:7)' * 1e9;
dps = tm_dt_to_dps(dt_1x6);

% covariance tensor
ct_1x21  = m(8:28)' * 1e18;
dps = tm_ct_to_dps(ct_1x21, dps);

% plot(h, S, 'k.');
% hold(h, 'on');
% plot(h, S_fit, 'r-');
% xlim(h, [0 1+numel(S_fit)]);
% ylim(h, [0 max(max(S_fit), max(S)) * 1.1]);
% title(axh, sprintf('MD = %1.2f, FA = %1.2f\nMK_I = %1.2f, MK_A = %1.2f\n  MEAS RANK = %i', ...
%     dps.MD, dps.FA, dps.MKi, dps.MKa, n));

title(axh, sprintf('MD = %1.2f, FA = %1.2f\nMK_I = %1.2f, MK_A = %1.2f\n ', ...
    dps.MD, dps.FA, dps.MKi, dps.MKa));
title_str = {...
    ['MD = ' num2str(dps.MD, 2) ' ',char(181),'m^2/ms   FA = ' num2str(dps.FA, 2) ]; ...
    ['\mu_2^{iso}/MD^2 = ' num2str(dps.MKi/3, 2) '   \mu_2^{aniso}/MD^2 = ' num2str(dps.MKa/3, 2)]};
    
title(axh, title_str)
axis(axh2,'off') 

if (0)
    
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
    
    surface(axh2,x,z,y,c);
    
    caxis(axh2, [0 1.5] );
    colormap(axh2, fliplr(jet(100)));
    axis(axh2, 'off')
    axis(axh2, 'tight');
    axis(axh2, 'equal');
    axis(axh2, 'vis3d');
    camlight(axh2, 0, 0);

    material([0.7 0.4 0.8 1]); 
    shading(axh2, 'interp');

end