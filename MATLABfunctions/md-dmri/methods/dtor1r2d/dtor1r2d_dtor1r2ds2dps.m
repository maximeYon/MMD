function dps = dtor1r2d_dtor1r2ds2dps(dps, dtor1r2ds, opt)
% function dps = dtor1r2d_dtor1r2ds2dps(dps, dtor1r2ds)

%Per-voxel statistical measures
nn = size(dtor1r2ds.w,4);

%DTD parameters
dps = dtd_dtds2dps(dps, dtor1r2ds);

%T1w and T2w signals
dps.s_te0010 = msf_notfinite2zero(sum(exp(-10e-3*dtor1r2ds.r2).*dtor1r2ds.w,4));
dps.s_te0020 = msf_notfinite2zero(sum(exp(-20e-3*dtor1r2ds.r2).*dtor1r2ds.w,4));
dps.s_te0050 = msf_notfinite2zero(sum(exp(-50e-3*dtor1r2ds.r2).*dtor1r2ds.w,4));
dps.s_te0100 = msf_notfinite2zero(sum(exp(-100e-3*dtor1r2ds.r2).*dtor1r2ds.w,4));
dps.s_tr0020 = msf_notfinite2zero(sum((1-exp(-20e-3*dtor1r2ds.r1)).*dtor1r2ds.w,4));
dps.s_tr0050 = msf_notfinite2zero(sum((1-exp(-50e-3*dtor1r2ds.r1)).*dtor1r2ds.w,4));
dps.s_tr0100 = msf_notfinite2zero(sum((1-exp(-100e-3*dtor1r2ds.r1)).*dtor1r2ds.w,4));
dps.s_tr0200 = msf_notfinite2zero(sum((1-exp(-200e-3*dtor1r2ds.r1)).*dtor1r2ds.w,4));
dps.s_tr0500 = msf_notfinite2zero(sum((1-exp(-500e-3*dtor1r2ds.r1)).*dtor1r2ds.w,4));
dps.s_tr1000 = msf_notfinite2zero(sum((1-exp(-1000e-3*dtor1r2ds.r1)).*dtor1r2ds.w,4));

% %Means
dps.mr1 = msf_notfinite2zero(sum(dtor1r2ds.r1.*dtor1r2ds.w,4)./dps.s0);
dps.mr2 = msf_notfinite2zero(sum(dtor1r2ds.r2.*dtor1r2ds.w,4)./dps.s0);

% %Variances
dps.vr1 = msf_notfinite2zero(sum((dtor1r2ds.r1-repmat(dps.mr1,[1 1 1 nn])).^2.*dtor1r2ds.w,4)./dps.s0);
dps.vr2 = msf_notfinite2zero(sum((dtor1r2ds.r2-repmat(dps.mr2,[1 1 1 nn])).^2.*dtor1r2ds.w,4)./dps.s0);

% %Covariances
dps.cvdisor1 = msf_notfinite2zero(sum((dtor1r2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtor1r2ds.r1-repmat(dps.mr1,[1 1 1 nn])).*dtor1r2ds.w,4)./dps.s0);
dps.cvdisor2 = msf_notfinite2zero(sum((dtor1r2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtor1r2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*dtor1r2ds.w,4)./dps.s0);
dps.cvsddeltar1 = msf_notfinite2zero(sum((dtor1r2ds.r1-repmat(dps.mr1,[1 1 1 nn])).*(dtor1r2ds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtor1r2ds.w,4)./dps.s0);
dps.cvsddeltar2 = msf_notfinite2zero(sum((dtor1r2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*(dtor1r2ds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtor1r2ds.w,4)./dps.s0);
dps.cvr1r2 = msf_notfinite2zero(sum((dtor1r2ds.r1-repmat(dps.mr1,[1 1 1 nn])).*(dtor1r2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*dtor1r2ds.w,4)./dps.s0);

% %Normalized measures
dps.nvr1 = msf_notfinite2zero(dps.vr1./dps.mr1.^2);
dps.nvr2 = msf_notfinite2zero(dps.vr2./dps.mr2.^2);
