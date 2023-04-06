function dps = dtod_dtods2dps(dps, dtods, opt)
% function dps = dtod_dtods2dps(dps, dtods)

%Per-voxel statistical measures
nn = size(dtods.w,4);

%DTD parameters
dps = dtd_dtds2dps(dps, dtods);

% %Means
% dps.mr1 = msf_notfinite2zero(sum(dtods.r1.*dtods.w,4)./dps.s0);
% dps.mr2 = msf_notfinite2zero(sum(dtods.r2.*dtods.w,4)./dps.s0);
% 
% %Variances
% dps.vr1 = msf_notfinite2zero(sum((dtods.r1-repmat(dps.mr1,[1 1 1 nn])).^2.*dtods.w,4)./dps.s0);
% dps.vr2 = msf_notfinite2zero(sum((dtods.r2-repmat(dps.mr2,[1 1 1 nn])).^2.*dtods.w,4)./dps.s0);
% 
% %Covariances
% dps.cvdisor1 = msf_notfinite2zero(sum((dtods.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtods.r1-repmat(dps.mr1,[1 1 1 nn])).*dtods.w,4)./dps.s0);
% dps.cvdisor2 = msf_notfinite2zero(sum((dtods.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtods.r2-repmat(dps.mr2,[1 1 1 nn])).*dtods.w,4)./dps.s0);
% dps.cvsddeltar1 = msf_notfinite2zero(sum((dtods.r1-repmat(dps.mr1,[1 1 1 nn])).*(dtods.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtods.w,4)./dps.s0);
% dps.cvsddeltar2 = msf_notfinite2zero(sum((dtods.r2-repmat(dps.mr2,[1 1 1 nn])).*(dtods.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtods.w,4)./dps.s0);
% dps.cvr1r2 = msf_notfinite2zero(sum((dtods.r1-repmat(dps.mr1,[1 1 1 nn])).*(dtods.r2-repmat(dps.mr2,[1 1 1 nn])).*dtods.w,4)./dps.s0);
% 
% %Normalized measures
% dps.nvr1 = msf_notfinite2zero(dps.vr1./dps.mr1.^2);
% dps.nvr2 = msf_notfinite2zero(dps.vr2./dps.mr2.^2);
