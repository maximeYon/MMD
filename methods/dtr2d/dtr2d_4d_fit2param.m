function dps = dtr2d_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtr2d_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
opt = dtr2d_opt(opt);
%dps = mdm_mfs_load(mfs_fn);

% Hack to allow mgui to access this function
if ischar(mfs_fn)
    dps = mdm_mfs_load(mfs_fn);
else
    m = mfs_fn;
    dps.m = m;
end


% create parameter maps and save them

m = dps.m;
dps = rmfield(dps,'m');
sz = size(m);

ind = false(sz(4),1);
ind(2:6:end) = 1;
nn = (sz(4)-1)/6;

dpar = m(:,:,:,circshift(ind,0,1));
dperp = m(:,:,:,circshift(ind,1,1));
theta = m(:,:,:,circshift(ind,2,1));
phi = m(:,:,:,circshift(ind,3,1));
r2 = m(:,:,:,circshift(ind,4,1));
w = m(:,:,:,circshift(ind,5,1));

%Calculate derived parameters
diso = (dpar + 2*dperp)/3;
daniso = (dpar - dperp)/3;
ddelta = msf_notfinite2zero(daniso./diso);
sdaniso = daniso.^2;
sddelta = msf_notfinite2zero(sdaniso./diso.^2);
dratio = msf_notfinite2zero(dpar./dperp);
[dxx,dyy,dzz,dxy,dxz,dyz] = dtr2d_pars2elements(dpar,dperp,theta,phi);

dtr2ds = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,'diso',diso,'daniso',daniso,'ddelta',ddelta,...
    'sdaniso',sdaniso,'sddelta',sddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,'dxy',dxy,'dxz',dxz,'dyz',dyz,'r2',r2);


    function dps = dtr2ds2dps(dps, dtr2ds)

        %Per-voxel statistical measures
        dps.s0 = sum(dtr2ds.w,4);

        %Means
        dps.mdiso = msf_notfinite2zero(sum(dtr2ds.diso.*dtr2ds.w,4)./dps.s0);
        dps.mdaniso = msf_notfinite2zero(sum(dtr2ds.daniso.*dtr2ds.w,4)./dps.s0);
        dps.msdaniso = msf_notfinite2zero(sum(dtr2ds.sdaniso.*dtr2ds.w,4)./dps.s0);
        dps.msddelta = msf_notfinite2zero(sum(dtr2ds.sddelta.*dtr2ds.w,4)./dps.s0);
        dps.mdxx = msf_notfinite2zero(sum(dtr2ds.dxx.*dtr2ds.w,4)./dps.s0);
        dps.mdyy = msf_notfinite2zero(sum(dtr2ds.dyy.*dtr2ds.w,4)./dps.s0);
        dps.mdzz = msf_notfinite2zero(sum(dtr2ds.dzz.*dtr2ds.w,4)./dps.s0);
        dps.mdxy = msf_notfinite2zero(sum(dtr2ds.dxy.*dtr2ds.w,4)./dps.s0);
        dps.mdxz = msf_notfinite2zero(sum(dtr2ds.dxz.*dtr2ds.w,4)./dps.s0);
        dps.mdyz = msf_notfinite2zero(sum(dtr2ds.dyz.*dtr2ds.w,4)./dps.s0);
        dps.mr2 = msf_notfinite2zero(sum(dtr2ds.r2.*dtr2ds.w,4)./dps.s0);

        %Variances
        dps.vdiso = msf_notfinite2zero(sum((dtr2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).^2.*dtr2ds.w,4)./dps.s0);
        dps.vsdaniso = msf_notfinite2zero(sum((dtr2ds.sdaniso-repmat(dps.msdaniso,[1 1 1 nn])).^2.*dtr2ds.w,4)./dps.s0);
        dps.vsddelta = msf_notfinite2zero(sum((dtr2ds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).^2.*dtr2ds.w,4)./dps.s0);
        dps.vr2 = msf_notfinite2zero(sum((dtr2ds.r2-repmat(dps.mr2,[1 1 1 nn])).^2.*dtr2ds.w,4)./dps.s0);

        %Covariances
        dps.cvdisosdaniso = msf_notfinite2zero(sum((dtr2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2ds.sdaniso-repmat(dps.msdaniso,[1 1 1 nn])).*dtr2ds.w,4)./dps.s0);
        dps.cvdisosddelta = msf_notfinite2zero(sum((dtr2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2ds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtr2ds.w,4)./dps.s0);
        dps.cvdisor2 = msf_notfinite2zero(sum((dtr2ds.diso-repmat(dps.mdiso,[1 1 1 nn])).*(dtr2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*dtr2ds.w,4)./dps.s0);
        dps.cvsddeltar2 = msf_notfinite2zero(sum((dtr2ds.r2-repmat(dps.mr2,[1 1 1 nn])).*(dtr2ds.sddelta-repmat(dps.msddelta,[1 1 1 nn])).*dtr2ds.w,4)./dps.s0);

        %Normalized measures
        dps.mdanison = msf_notfinite2zero(dps.mdaniso./dps.mdiso);
        dps.msdanison = msf_notfinite2zero(dps.msdaniso./dps.mdiso.^2);
        dps.vdison = msf_notfinite2zero(dps.vdiso./dps.mdiso.^2);
        dps.vsdanison = msf_notfinite2zero(dps.vsdaniso./dps.msdaniso.^2);
        dps.vsddeltan = msf_notfinite2zero(dps.vsddelta./dps.msddelta.^2);
        dps.vr2n = msf_notfinite2zero(dps.vr2./dps.mr2.^2);


        dps.Vl = 5/2 * 4/5*dps.msdaniso;
        dps.MKi = 3 * dps.vdison; % Multiply by 3 to get kurtosis
        dps.MKa = 3 * 4/5*dps.msdanison;

        % Calculate uFA. Take real component to avoid complex values due to
        % sqrt of negative variances.
        dps.ufa_old = real(sqrt(3/2) * sqrt(1./(dps.mdiso.^2./dps.Vl+1))); % Lasic (2014)
        dps.ufa     = real(sqrt(3/2) * sqrt( dps.Vl ./ (dps.Vl + dps.vdiso + dps.mdiso.^2) )); % Szczepankiewicz (2016)
        
        dps.signaniso = sign(dps.mdanison);
    end

dps = dtr2ds2dps(dps, dtr2ds);

%Per-bin statistical measures
for nbin = 1:numel(opt.dtr2d.bin_disomax)
    ind_bin = false([sz(1) sz(2) sz(3) nn 4]);
    ind_bin(:,:,:,:,1) = diso >= opt.dtr2d.bin_disomin(nbin);
    ind_bin(:,:,:,:,2) = diso <= opt.dtr2d.bin_disomax(nbin);
    ind_bin(:,:,:,:,3) = dratio >= opt.dtr2d.bin_dratiomin(nbin);
    ind_bin(:,:,:,:,4) = dratio <= opt.dtr2d.bin_dratiomax(nbin);
    ind_bin = all(ind_bin,5);

    dps_bin.no = nbin;
    dtr2ds_temp = dtr2ds;
    dtr2ds_temp.w = dtr2ds.w.*ind_bin;
    dps_bin = dtr2ds2dps(dps_bin, dtr2ds_temp);
    dps.bin{nbin} = dps_bin;
end

if (~isempty(dps_fn)) mdm_dps_save(dps, dps.s, dps_fn, opt); end

end

