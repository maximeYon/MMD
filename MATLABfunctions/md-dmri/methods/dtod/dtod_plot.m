function m = dtod_plot(S, xps, axh, axh2, opt)
% function dtod_plot(S, xps, axh, axh2)

if (nargin < 4), axh2 = []; end
if (nargin < 5), opt = []; end

method = 'dtod';

opt = mdm_opt(opt);
opt = dtd_opt(opt);
opt = feval([method '_opt'],opt);
opt = mplot_opt(opt);

opt.(method).dmin = .01/max(xps.b);
opt.(method).dmax = 1/min(xps.b(xps.b>1e7));
opt.(method).rmin = .01*min(xps.momega(xps.momega>0)/2/pi);
opt.(method).rmax = 100*max(xps.momega/2/pi);
opt.(method).n = 10;
opt.(method).n_proliferation = 20;
opt.(method).n_darwin = 20;
opt.(method).maps_omega = mean(xps.momega);
%     opt.(method)
    
S = abs(S);

% Show signal and fit
m = mplot_signal_and_fit(S, xps, @dtod_1d_data2fit, @dtod_1d_fit2data, axh, opt);
title(axh,['<\omega> / 2\pi = ' num2str(min(xps.momega(xps.b>.5*max(xps.b))/2/pi),3) ' - ' num2str(max(xps.momega(xps.b>.5*max(xps.b))/2/pi),3) ' Hz'])

% Get dps from dtod_4d_fit2param
mfs.m = zeros(1,1,1,numel(m)); mfs.m(1,1,1,:) = m;
dps = dtod_4d_fit2param(mfs.m,[],opt);

[dpar,dperp,theta,phi,w] = dtod_4d_m2parso(mfs.m,mean(xps.momega));
[dxx,dyy,dzz,dxy,dxz,dyz] = dtd_pars2elements(squeeze(dpar),squeeze(dperp),squeeze(theta),squeeze(phi));
dtd_1x6 = [dxx dyy dzz sqrt(2)*[dxy dxz dyz]];

% Plot the tensor distribution
% [dtd_1x6,w] = dtd_dist2nx6w(dtd_m2dtd(m));
mplot_dtd(dtd_1x6, w, opt.dtod.dmin, opt.dtod.dmax, axh2, opt);
mplot_dtd_addstats(dps, axh2, opt);
mplot_dtd_addtitle(dps, axh2, opt);

