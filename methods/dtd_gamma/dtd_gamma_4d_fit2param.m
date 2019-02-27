function dps = dtd_gamma_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_gamma_4d_fit2param(mfs_fn, dps_fn, opt)
%
% Calculates derived parameters from the primary parameters of the
% gamma fit.
%
% Parameters introduced in Lasic et al, Front. Phys. 2, 11 (2014). http://dx.doi.org/10.3389/fphy.2014.00011.
% Renamed in Westin et al, Neuroimage 135, 345 (2016), Szczepankiewicz et
% al, Neuroimage 142, 522 (2016), and Topgaard. J. Magn. Reson. 275, 98 (2017).

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
%dps = mdm_mfs_load(mfs_fn);

% Hack to allow mgui to access this function
if ischar(mfs_fn)
    dps = mdm_mfs_load(mfs_fn);
else
    m = mfs_fn;
    dps.m = m;
end

dps.s0 = dps.m(:,:,:,1);

% Parameters according to terminology in Szczepankiewicz (2016)
dps.MD = dps.m(:,:,:,2)*1e9; % mean diffusivity
dps.Vi = dps.m(:,:,:,3)*1e9*1e9; % \mu_2^{iso} in Lasic (2014)
dps.Va = dps.m(:,:,:,4)*1e9*1e9; % \Delta\mu_2 in Lasic (2014)
dps.Vt = dps.Vi + dps.Va; % \mu_2 in Lasic (2014)

dps.NVi = dps.Vi ./ dps.MD.^2; % \tilde{\mu}_2^{iso} in Lasic (2014), C_{MD} in Westin (2016)
dps.NVa = dps.Va ./ dps.MD.^2; % \Delta\tilde{\mu}_2 in Lasic (2014)
dps.NVt = dps.Vt ./ dps.MD.^2; % \tilde{\mu}_2 in Lasic (2014)

dps.MKi = 3 * dps.NVi; % Multiply by 3 to get kurtosis
dps.MKa = 3 * dps.NVa;
dps.MKt = 3 * dps.NVt;

dps.Vl = 5/2 * dps.Va;

% Topgaard. J. Magn. Reson. 275, 98 (2017). https://dx.doi.org/10.1016/j.jmr.2016.12.007
% Recommened for comparison with results from the dtd method.
dps.mdiso = dps.m(:,:,:,2); % Mean isotropic diffusivity, see Eq. (68)
dps.mu2iso = dps.m(:,:,:,3);
dps.mu2aniso = dps.m(:,:,:,4);
dps.mu2ison = dps.mu2iso./dps.mdiso.^2; % Normalized
dps.mu2anison = dps.mu2aniso./dps.mdiso.^2;

dps.vdiso = dps.mu2iso;  % Variance of isotropic diffusivities, see Eqs. (48) and (69)
dps.msdaniso = 5/4*dps.mu2aniso;  % Mean-square anisotropic diffusivity, see Eqs. (50) and (69)
dps.vdison = dps.vdiso./dps.mdiso.^2; % Normalized
dps.msdanison = dps.msdaniso./dps.mdiso.^2;

% Calculate uFA. Take real component to avoid complex values due to
% sqrt of negative variances.
dps.ufa_old = real(sqrt(3/2) * sqrt(1./(dps.MD.^2./dps.Vl+1))); % Lasic (2014)
dps.ufa     = real(sqrt(3/2) * sqrt( dps.Vl ./ (dps.Vl + dps.Vi + dps.MD.^2) )); % Szczepankiewicz (2016)
dps.uFA = dps.ufa;

for i = 5:size(dps.m, 4)
    nam = ['s' num2str(i-4)];
    dps.(nam) = dps.m(:,:,:,i);
end


if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end
