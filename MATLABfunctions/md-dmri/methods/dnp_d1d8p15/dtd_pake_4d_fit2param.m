function dps = dtd_pake_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_pake_4d_fit2param(dps_fn, dps_fn, opt)

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

% create parameter maps and save them
dps.s0      = dps.m(:,:,:,1);
dps.diso     = dps.m(:,:,:,2);
dps.ddelta   = dps.m(:,:,:,3);

dps.dpar             = dps.diso.*(1 + 2*dps.ddelta);
dps.dperp            = dps.diso.*(1 - dps.ddelta);

dps.vlambda         = 2*(dps.diso.*dps.ddelta).^2;
dps.ufa             = sqrt(3/2)*sqrt(1./(dps.diso.^2./dps.vlambda+1));
dps.cmu             = dps.ufa.^2;

dps.mdiso = dps.diso;  % Mean-square anisotropic diffusivity, see Eqs. (50) and (69)
dps.vdiso = zeros(size(dps.s0));  % Variance of isotropic diffusivities, see Eqs. (48) and (69)
dps.msdaniso = (dps.diso.*dps.ddelta).^2;  % Mean-square anisotropic diffusivity, see Eqs. (50) and (69)
dps.nvdiso = dps.vdiso./dps.mdiso.^2; % Normalized
dps.nmsdaniso = dps.msdaniso./dps.mdiso.^2;

dps.MKi = 3 * dps.nvdiso; % Multiply by 3 to get kurtosis
dps.MKa = 3 * 4/5*dps.nmsdaniso;

dps.signaniso = sign(dps.ddelta);

if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end



