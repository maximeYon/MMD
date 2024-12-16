function dps = dtr1r2d_dtr1r2ds2dps(dps, dtr1r2ds)
% function dps = dtr1r2d_dtr1r2ds2dps(dps, dtr1r2ds)

%Per-voxel statistical measures
nn = size(dtr1r2ds.w,4);

%DTD parameters
dps = dtd_dtds2dps(dps, dtr1r2ds);

%Means
dps.mr1 = msf_notfinite2zero(sum(dtr1r2ds.r1.*dtr1r2ds.w,4)./dps.s0);
dps.mr2 = msf_notfinite2zero(sum(dtr1r2ds.r2.*dtr1r2ds.w,4)./dps.s0);

%Variances
dps.vr1 = msf_notfinite2zero(sum((dtr1r2ds.r1-repmat(dps.mr1,[1 1 1 nn])).^2.*dtr1r2ds.w,4)./dps.s0);
dps.vr2 = msf_notfinite2zero(sum((dtr1r2ds.r2-repmat(dps.mr2,[1 1 1 nn])).^2.*dtr1r2ds.w,4)./dps.s0);

%Covariances
dps.cvdisor1 = msf_notfinite2zero(sum((dtr1r2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr1r2ds.r1-repmat(dps.mr1,[1 1 1 nn])).*dtr1r2ds.w,4)./dps.s0);
dps.cvdisor2 = msf_notfinite2zero(sum((dtr1r2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr1r2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*dtr1r2ds.w,4)./dps.s0);
dps.cvsddeltar1 = msf_notfinite2zero(sum((dtr1r2ds.r1-repmat(dps.mr1,[1 1 1 nn])).*(dtr1r2ds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtr1r2ds.w,4)./dps.s0);
dps.cvsddeltar2 = msf_notfinite2zero(sum((dtr1r2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*(dtr1r2ds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtr1r2ds.w,4)./dps.s0);
dps.cvr1r2 = msf_notfinite2zero(sum((dtr1r2ds.r1-repmat(dps.mr1,[1 1 1 nn])).*(dtr1r2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*dtr1r2ds.w,4)./dps.s0);

%Normalized measures
dps.nvr1 = msf_notfinite2zero(dps.vr1./dps.mr1.^2);
dps.nvr2 = msf_notfinite2zero(dps.vr2./dps.mr2.^2);

%T1w and T2w signals
dps.s_te0010 = msf_notfinite2zero(sum(exp(-10e-3*dtr1r2ds.r2).*dtr1r2ds.w,4));
dps.s_te0020 = msf_notfinite2zero(sum(exp(-20e-3*dtr1r2ds.r2).*dtr1r2ds.w,4));
dps.s_te0050 = msf_notfinite2zero(sum(exp(-50e-3*dtr1r2ds.r2).*dtr1r2ds.w,4));
dps.s_te0100 = msf_notfinite2zero(sum(exp(-100e-3*dtr1r2ds.r2).*dtr1r2ds.w,4));
dps.s_tr0020 = msf_notfinite2zero(sum((1-exp(-20e-3*dtr1r2ds.r1)).*dtr1r2ds.w,4));
dps.s_tr0050 = msf_notfinite2zero(sum((1-exp(-50e-3*dtr1r2ds.r1)).*dtr1r2ds.w,4));
dps.s_tr0100 = msf_notfinite2zero(sum((1-exp(-100e-3*dtr1r2ds.r1)).*dtr1r2ds.w,4));
dps.s_tr0200 = msf_notfinite2zero(sum((1-exp(-200e-3*dtr1r2ds.r1)).*dtr1r2ds.w,4));
dps.s_tr0500 = msf_notfinite2zero(sum((1-exp(-500e-3*dtr1r2ds.r1)).*dtr1r2ds.w,4));
dps.s_tr1000 = msf_notfinite2zero(sum((1-exp(-1000e-3*dtr1r2ds.r1)).*dtr1r2ds.w,4));
dps.s_tr2000 = msf_notfinite2zero(sum((1-exp(-2000e-3*dtr1r2ds.r1)).*dtr1r2ds.w,4));
dps.s_tr5000 = msf_notfinite2zero(sum((1-exp(-5000e-3*dtr1r2ds.r1)).*dtr1r2ds.w,4));

%T1w and T2w means
[tr,te] = ndgrid([1e3-1e-3,4,2,1],[0,25,50,100]*1e-3);
tr = tr(:); te = te(:);
for ntr = 1:numel(tr)
    w_temp = (1-exp(-tr(ntr)*dtr1r2ds.r1)).*exp(-te(ntr)*dtr1r2ds.r2).*dtr1r2ds.w;
    trte_str = strcat('_tr',sprintf('%06d',tr(ntr)*1e3),'_te',sprintf('%06d',te(ntr)*1e3));
    dps.(strcat('mdiso',trte_str)) = msf_notfinite2zero(sum(dtr1r2ds.diso.*w_temp,4)./sum(w_temp,4));
    dps.(strcat('msddelta',trte_str)) = msf_notfinite2zero(sum(dtr1r2ds.sddelta.*w_temp,4)./sum(w_temp,4));
end

