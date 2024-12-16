function dps = dtr2d_dtr2ds2dps(dps, dtr2ds)
% function dps = dtr2d_dtr2ds2dps(dps, dtr2ds)

%Per-voxel statistical measures
nn = size(dtr2ds.w,4);

%DTD parameters
dps = dtd_dtds2dps(dps, dtr2ds);

%T2-weighted signals
dps.s_te10 = msf_notfinite2zero(sum(exp(-10e-3*dtr2ds.r2).*dtr2ds.w,4));
dps.s_te20 = msf_notfinite2zero(sum(exp(-20e-3*dtr2ds.r2).*dtr2ds.w,4));
dps.s_te50 = msf_notfinite2zero(sum(exp(-50e-3*dtr2ds.r2).*dtr2ds.w,4));
dps.s_te100 = msf_notfinite2zero(sum(exp(-100e-3*dtr2ds.r2).*dtr2ds.w,4));

%ADC vs. b and TE
b = .2e9; TE = 0e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b0200_te000 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = .5e9; TE = 0e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b0500_te000 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 0e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te000 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 1e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te001 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 10e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te010 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 20e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te020 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 25e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te025 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 30e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te030 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 35e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te035 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 40e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te040 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 45e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te045 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 50e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te050 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;
b = 1e9; TE = 100e-3;
s0 = msf_notfinite2zero(sum(exp(-0e9*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
sw = msf_notfinite2zero(sum(exp(-b*dtr2ds.diso - TE*dtr2ds.r2).*dtr2ds.w,4));
dps.adc_b1000_te100 = msf_notfinite2zero((log(s0) - log(sw))./b)/1e-9;

%Means
dps.mr2 = msf_notfinite2zero(sum(dtr2ds.r2.*dtr2ds.w,4)./dps.s0);

%Variances
dps.vr2 = msf_notfinite2zero(sum((dtr2ds.r2-repmat(dps.mr2,[1 1 1 nn])).^2.*dtr2ds.w,4)./dps.s0);

%Covariances
dps.cvdisor2 = msf_notfinite2zero(sum((dtr2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*dtr2ds.w,4)./dps.s0);
dps.cvsddeltar2 = msf_notfinite2zero(sum((dtr2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*(dtr2ds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtr2ds.w,4)./dps.s0);

%Normalized measures
dps.nvr2 = msf_notfinite2zero(dps.vr2./dps.mr2.^2);
