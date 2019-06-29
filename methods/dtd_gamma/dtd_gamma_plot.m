function dtd_gamma_plot(S, xps, axh, axh2, opt)
% function dtd_gamma_plot(S, xps, h, h2)

if (nargin < 3), axh  = gca; end
if (nargin < 4), axh2 = []; end
if (nargin < 5), opt = []; end

opt = mdm_opt(opt);
opt = dtd_gamma_opt(opt);
opt = dtd_opt(opt);
opt = mplot_opt(opt);

if any(~isreal(S(:)))
    S = abs(S);
end

m = mplot_s_vs_b_by_b_delta(S, xps, ...
    @dtd_gamma_1d_data2fit, @dtd_gamma_1d_fit2data, @dtd_gamma_opt, axh);

% Get dps from dtd_gamma_4d_fit2param
mfs.m = zeros(1,1,1,numel(m)); mfs.m(1,1,1,:) = m;
dps = dtd_gamma_4d_fit2param(mfs.m);

% Plot the tensor distribution
[dtd_1x6,w] = dtd_dist2nx6w(dtd_m2dtd([1 1e-13 1e-13 0 0 1])); %dummy distribution to get axes
mplot_dtd(dtd_1x6, w, opt.dtd.dmin, opt.dtd.dmax, axh2, opt);
mplot_dtd_addstats(dps, axh2, opt);
mplot_dtd_addtitle(dps, axh2, opt);
