function dps = tm_ct_to_dps(ct, dps, f_reshape, reg)
% function dps = tm_ct_to_dps(dt, dps, f_reshape)
%
% ct  - diffusion tensor covariance 4th order tensor(s) in 1x21 format
% dps - display parameter structure
%
% optional
% f_reshape - function to reshape output parameters
% reg       - regularization parameter
%
% output:
% dps.(V_MD, V_iso, VM_MD1, V_iso1, V_shear, V_shear1, ...
%        C_MD, C_mu, C_M, C_c, K_bulk, K_shear, MK, K_mu)

if (nargin < 2), dps = []; end
if (nargin < 3) || (isempty(f_reshape)), f_reshape = @(x,n) x; end
if (nargin < 4), reg = eps; end

% fourth order tensor operators needed for the covariance calculations
[E_bulk, E_shear, E_iso] = tm_1x21_iso();

% DTD paramters derived from first and second cumulant term
% ********* I.V_MD = I.V_MD1 - I.V_MD2 *************
dps.V_MD   = f_reshape(tm_inner(ct, E_bulk), 1);       %< C ,   E_bulk>
dps.V_iso  = f_reshape(tm_inner(ct, E_iso), 1);        %< C ,   E_iso>

dps.V_MD1  = dps.V_MD  + dps.V_MD2; % V_MD2 comes from tm_dt_to_dps
dps.V_iso1 = dps.V_iso + dps.V_iso2;

% ***** I.V_shear = I.V_shear1 - I.V_shear2 *******
dps.V_shear  = f_reshape(tm_inner(ct, E_shear), 1);     %< C ,    E_shear>
dps.V_shear1 = dps.V_shear + dps.V_shear2; % V_shear2 comes from tm_dt_to_dps

% ********* Normalized variance measures ************
dps.C_MD    = dps.V_MD ./ max(dps.V_MD1,reg);
dps.C_mu    = 1.5 * dps.V_shear1 ./ max(dps.V_iso1, reg);
dps.C_M     = 1.5 * dps.V_shear2 ./ max(dps.V_iso2, reg);
dps.C_c     = dps.C_M ./ max(dps.C_mu, reg);

% ********* Kurtosis measures ************
% Naming these according to the dtd_gamma nomenclature
dps.MKi  = 3 * dps.V_MD ./ max(dps.V_MD2, reg);           % 
dps.MKa  = (6/5) * dps.V_shear1 ./ max(dps.V_MD2, reg);   % K_micro in Westin16
dps.MKt  = dps.MKi + dps.MKa;                   % 
dps.MKad = (6/5) * dps.V_shear ./ max(dps.V_MD2, reg);    % anisotropy and dispersion
dps.MK   = dps.MKad + dps.MKi;                  % conventional kurtosis
dps.MKd  = dps.MKa - dps.MKad;                  % conventional kurtosis
dps.uFA  = sqrt(dps.C_mu);


dps.S_I = sqrt(dps.V_MD .* (dps.V_MD > 0));
dps.S_A = sqrt(dps.V_shear1 .* (dps.V_shear1 > 0));




