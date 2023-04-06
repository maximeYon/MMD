% Init
% rng(1);

% Generate a diffusion tensor distribution
ad_1 = 1.7e-9; % axial diffusivity
rd_1 = 0.1e-9; % radial diffusivity
w_1  = 0.5;

ad_2 = 2.0e-9;
rd_2 = 1.0e-9;
w_2 = 0.3;

ad_3 = 3.0e-9;
rd_3 = ad_3;
w_3  = 0.1;

n_all = uvec_elstat(512); % tensors point all across the sphere

dt = cat(1, ...
    tm_1x3_to_1x6(ad_1, rd_1, n_all), ...
    tm_1x3_to_1x6(ad_2, rd_2, n_all), ...
    tm_1x3_to_1x6(ad_3, rd_3, n_all));

f = @(x) repmat(x, size(n_all, 1), 1);
w = cat(1, f(w_1), f(w_2), f(w_2));
w = w / sum(w);

% for plotting
n = [1 0 0];

dt_native = cat(1, ...
    tm_1x3_to_1x6(ad_1, rd_1, n), ...
    tm_1x3_to_1x6(ad_2, rd_2, n), ...
    tm_1x3_to_1x6(ad_3, rd_3, n));

w_native = [w_1 w_2 w_3]; w = w / sum(w);

% Prepare the measurement
b_max = 8e9;
b = linspace(0,1, 8).^2 * b_max * 0.99 + 0.01 * b_max;

b_delta = [1 0 -0.5]; % tensor shapes

xps = {};
for c_delta = 1:numel(b_delta)
    for c_b = 1:numel(b)
        n_rand = randi(9) + 20;
        
        n = uvec_tricosa;
        
        ab = b(c_b) * (b_delta(c_delta) + 0.5) / 1.5;
        rb = (b(c_b) - ab) / 2;
        
        bt = tm_1x3_to_1x6(ab, rb, n); % b-tensor
        
        % prepare eXperiment Parameter Structure
        xps{end+1} = mdm_xps_from_bt(bt); 
        
    end
end
xps = mdm_xps_merge(xps);

% Generate data
S = exp(-xps.bt * dt') * w + randn(xps.n,1) * 0.001;

% Plot
clf; 
axes('position', [0.15 0.17 0.8 0.75]);

opt = dtd_pa_opt();
opt.mplot.fs = 14;

dmin = 10^(-10.2);
dmax = 10^(-8.2);

opt.mplot.dtd_col_mode = 0;

c_case = 2;
switch (c_case)
    case 1 % Underlying DTD
        mplot_dtd(dt_native, w_native, dmin, dmax, gca, opt);
    
    case 2 % Signal and fit
        mplot_s_vs_b_by_b_delta(S, xps, @dtd_pa_1d_data2fit, @dtd_pa_1d_fit2data, opt, gca, []);
        ylim([0.01 1]);

    case 3 % Fitted distribution
        m = dtd_pa_1d_data2fit(S, xps, dtd_pa_opt());
        [dtd_1x6, w] = dtd_pa_dtd_to_t1x6(dtd_pa_m2dtd(m));
        mplot_dtd(dtd_1x6, w, dmin, dmax, gca, opt);
end



