function dps = dtd_dtd2dps(dtd)
% dps = dtd_dtd2dps(dtd)
%
% Converts discrete diffusion tensor distribution to derived tensors and
% parameters. Equation numbers refer to definitions given in
% Topgaard. J. Magn. Reson. 275, 98 (2017).
% http://dx.doi.org/10.1016/j.jmr.2016.12.007
%
% Input arguments
% dtd: diffusion tensor distribution
%
% Output
% dps: derived parameters structure with fields
% s0: total amplitude = b0 signal
% miso: mean isotropic diffusivity, Eq (7.26)
% viso: variance of isotropic diffusivities, Eq (7.27)
% maniso: mean anisotropic diffusivity, Eq (7.28)
% vaniso: variance of anistropic diffusivities
% msqaniso: mean-square anisotropic diffusivity, Eq (7.29)
% vsqaniso: variance of squared anisotropic diffusivities
% t1x6: mean diffusion tensor in Voigt notation, 
% ct1x21: covariance tensor in Voigt-like notation, Eq (7.78)

dps.s0 = zeros([1 1]);
dps.t1x6 = zeros([1 6]);
dps.lambdazzvec = zeros([1 3]);
dps.lambdaxxvec = zeros([1 3]);
dps.lambdayyvec = zeros([1 3]);
dps.lambda11vec = zeros([1 3]);
dps.lambda22vec = zeros([1 3]);
dps.lambda33vec = zeros([1 3]);
dtiparam = {'trace','iso','lambda33','lambda22','lambda11','lambdazz','lambdaxx','lambdayy','vlambda',...
    'delta','eta','s','p','l','fa','cs','cl','cp','cm'};
param = {dtiparam{:},'miso','viso','maniso','vaniso','msqaniso','vsqaniso'};
for nparam = 1:numel(param)
    eval(['dps.' param{nparam} ' = zeros([1 1 1]);']);
end

[n,dpar,dperp,theta,phi,w] = dtd_dist2par(dtd);
%Calculate derived parameters
diso = (dpar + 2*dperp)/3;
daniso = (dpar - dperp)/3;
ddelta = msf_notfinite2zero(daniso./diso);
sdaniso = daniso.^2;
sddelta = msf_notfinite2zero(sdaniso./diso.^2);
dratio = msf_notfinite2zero(dpar./dperp);

if n > 0
    %Per-voxel statistical measures
    %Means
    s0 = ones(1,n)*w;
    mdiso = diso'*w/s0;
    msdaniso = sdaniso'*w/s0;
    msddelta = sddelta'*w/s0;

    vdiso = (diso-mdiso)'.^2*w/s0;
    vsdaniso = (sdaniso-msdaniso)'.^2*w/s0; 
    vsddelta = (sddelta-msddelta)'.^2*w/s0;

        dps.mdiso = msf_notfinite2zero(sum(dtds.diso.*dtds.w,4)./dps.s0);
        dps.msdaniso = msf_notfinite2zero(sum(dtds.sdaniso.*dtds.w,4)./dps.s0);
        dps.msddelta = msf_notfinite2zero(sum(dtds.sddelta.*dtds.w,4)./dps.s0);
        dps.mdxx = msf_notfinite2zero(sum(dtds.dxx.*dtds.w,4)./dps.s0);
        dps.mdyy = msf_notfinite2zero(sum(dtds.dyy.*dtds.w,4)./dps.s0);
        dps.mdzz = msf_notfinite2zero(sum(dtds.dzz.*dtds.w,4)./dps.s0);
        dps.mdxy = msf_notfinite2zero(sum(dtds.dxy.*dtds.w,4)./dps.s0);
        dps.mdxz = msf_notfinite2zero(sum(dtds.dxz.*dtds.w,4)./dps.s0);
        dps.mdyz = msf_notfinite2zero(sum(dtds.dyz.*dtds.w,4)./dps.s0);
        dps.t1x6 = [dps.mdxx dps.mdyy dps.mdzz srt(2)*[dps.mdxy dps.mdxz dps.mdyz]];

        %Variances
        dps.vdiso = msf_notfinite2zero(sum((dtds.diso-repmat(dps.mdiso,[1 1 1 nn])).^2.*dtds.w,4)./dps.s0);
        dps.vsdaniso = msf_notfinite2zero(sum((dtds.sdaniso-repmat(dps.msdaniso,[1 1 1 nn])).^2.*dtds.w,4)./dps.s0);
        dps.vsddelta = msf_notfinite2zero(sum((dtds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).^2.*dtds.w,4)./dps.s0);

        %Covariances
        dps.cvdisosdaniso = msf_notfinite2zero(sum((dtds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtds.sdaniso-repmat(dps.msdaniso,[1 1 1 nn])).*dtds.w,4)./dps.s0);
        dps.cvdisosddelta = msf_notfinite2zero(sum((dtds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtds.w,4)./dps.s0);

        %Normalized measures
        dps.vdison = msf_notfinite2zero(dps.vdiso./dps.mdiso.^2);
        dps.vsdanison = msf_notfinite2zero(dps.vsdaniso./dps.msdaniso.^2);
        dps.vsddeltan = msf_notfinite2zero(dps.vsddelta./dps.msddelta.^2);
    
    %Normalized with miso
    viso_n = viso./miso^2;
    msqaniso_n = msqaniso./miso^2;

    dps.s0 = s0;
    dps.miso = miso;
    dps.viso = viso;
    dps.maniso = maniso;
    dps.vaniso = vaniso;
    dps.msqaniso = msqaniso;
    dps.vsqaniso = vsqaniso;

    [dtd_nx6,w] = dtd_dist2nx6w(dtd);
    dt1x6 = (dtd_nx6'*w)'/s0;
    dt3x3 = tm_1x6_to_3x3(dt1x6);
    dt = tm_3x3_to_tpars(dt3x3);

    dps.t1x6(1,:) = dt.t1x6;
    dps.lambdazzvec(1,:) = dt.lambdazzvec;
    dps.lambdaxxvec(1,:) = dt.lambdaxxvec;
    dps.lambdayyvec(1,:) = dt.lambdayyvec;
    dps.lambda11vec(1,:) = dt.lambda11vec;
    dps.lambda22vec(1,:) = dt.lambda22vec;
    dps.lambda33vec(1,:) = dt.lambda33vec;
    for nparam = 1:numel(dtiparam)
        eval(['dps.' dtiparam{nparam} '(1) = dt.' dtiparam{nparam} ';']);
    end

    dt1x21 = tm_1x6_to_1x21(dt1x6);
    dtd_nx21 = tm_1x6_to_1x21(dtd_nx6);
    dps.ct1x21 = (dtd_nx21'*w)'/s0 - dt1x21;
    
end
    