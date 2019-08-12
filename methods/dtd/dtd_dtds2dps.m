function dps = dtd_dtds2dps(dps, dtds)
% function dps = dtd_dtds2dps(dps, dtds)

%Per-voxel statistical measures
nn = size(dtds.w,4);

dps.s0 = sum(dtds.w,4);

%Means
dps.mdiso = msf_notfinite2zero(sum(dtds.diso.*dtds.w,4)./dps.s0);
dps.mdaniso = msf_notfinite2zero(sum(dtds.daniso.*dtds.w,4)./dps.s0);
dps.msdaniso = msf_notfinite2zero(sum(dtds.sdaniso.*dtds.w,4)./dps.s0);
dps.msddelta = msf_notfinite2zero(sum(dtds.sddelta.*dtds.w,4)./dps.s0);
dps.mdxx = msf_notfinite2zero(sum(dtds.dxx.*dtds.w,4)./dps.s0);
dps.mdyy = msf_notfinite2zero(sum(dtds.dyy.*dtds.w,4)./dps.s0);
dps.mdzz = msf_notfinite2zero(sum(dtds.dzz.*dtds.w,4)./dps.s0);
dps.mdxy = msf_notfinite2zero(sum(dtds.dxy.*dtds.w,4)./dps.s0);
dps.mdxz = msf_notfinite2zero(sum(dtds.dxz.*dtds.w,4)./dps.s0);
dps.mdyz = msf_notfinite2zero(sum(dtds.dyz.*dtds.w,4)./dps.s0);

%Variances
dps.vdiso = msf_notfinite2zero(sum((dtds.diso-repmat(dps.mdiso,[1 1 1 nn])).^2.*dtds.w,4)./dps.s0);
dps.vsdaniso = msf_notfinite2zero(sum((dtds.sdaniso-repmat(dps.msdaniso,[1 1 1 nn])).^2.*dtds.w,4)./dps.s0);
dps.vsddelta = msf_notfinite2zero(sum((dtds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).^2.*dtds.w,4)./dps.s0);

%Covariances
dps.cvdisosdaniso = msf_notfinite2zero(sum((dtds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtds.sdaniso-repmat(dps.msdaniso,[1 1 1 nn])).*dtds.w,4)./dps.s0);
dps.cvdisosddelta = msf_notfinite2zero(sum((dtds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtds.w,4)./dps.s0);

%Normalized measures
dps.nmdaniso = msf_notfinite2zero(dps.mdaniso./dps.mdiso);
dps.nmsdaniso = msf_notfinite2zero(dps.msdaniso./dps.mdiso.^2);
dps.nvdiso = msf_notfinite2zero(dps.vdiso./dps.mdiso.^2);
dps.nvsdaniso = msf_notfinite2zero(dps.vsdaniso./dps.msdaniso.^2);
dps.nvsddelta = msf_notfinite2zero(dps.vsddelta./dps.msddelta.^2);


dps.Vl = 5/2 * 4/5*dps.msdaniso;
dps.MKi = 3 * dps.nvdiso; % Multiply by 3 to get kurtosis
dps.MKa = 3 * 4/5*dps.nmsdaniso;

% Calculate uFA. Take real component to avoid complex values due to
% sqrt of negative variances.
dps.ufa_old = real(sqrt(3/2) * sqrt(1./(dps.mdiso.^2./dps.Vl+1))); % Lasic (2014)
dps.ufa     = real(sqrt(3/2) * sqrt( dps.Vl ./ (dps.Vl + dps.vdiso + dps.mdiso.^2) )); % Szczepankiewicz (2016)

dps.signaniso = sign(dps.nmdaniso);
