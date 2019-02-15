function dtd_plot(S, xps, axh, axh2)
% function dtd_plot(S, xps, axh, axh2)

if (nargin < 4), axh2 = []; end

opt = mdm_opt();
opt = dtd_opt(opt);
opt = mplot_opt(opt);

% Customize options
opt.dtd.dmin = .02/max(xps.b);
%[opt.dtd.dmin opt.dtd.dmax]

S = abs(S);

% Show signal and fit
m = mplot_signal_and_fit(S, xps, @dtd_1d_data2fit, @dtd_1d_fit2data, axh, opt);

% Plot the tensor distribution
[dtd_1x6,w] = dtd_dist2nx6w(dtd_m2dtd(m));
mplot_dtd(dtd_1x6, w, opt.dtd.dmin, opt.dtd.dmax, axh2, opt);
    
dtd = dtd_m2dtd(m);
dps1d = dtd_dist2dps1d(dtd);