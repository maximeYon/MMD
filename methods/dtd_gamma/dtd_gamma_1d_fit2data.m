function s = dtd_gamma_1d_fit2data(m, xps)
% function s = dtd_gamma_1d_fit2data(m, xps)

% convert to readable parameters
s0        = m(1);         % Signal at no diff encoding (baseline signal)
d_iso     = m(2);         % Average isotropic diffusivity (mean diffusivity)
mu2_iso   = m(3);         % Isotropic variance
mu2_aniso = m(4);         % Anisotropic variance

rs        = [1 m(5:end)]; % Relative signal (Series 1, 2, 3, ...)

% signal baseline coefficient that adjusts for variable baseline signal
% across multiple series (xps.s_ind).
sw = s0 * sum(  mtimes(  ones(size(xps.b)), rs    ) .* ...
                      (  xps.s_ind == 1:numel(rs) )  , 2);
                               
% Total diffusional variance (mu2) is the weighted sum of the isotropic 
% (mu2_iso) and anisotorpic (mu2_aniso) components.
mu2 = mu2_iso + mu2_aniso.*xps.b_delta.^2.*(xps.b_eta.^2+3)/3;

% Singnal equation based on the laplace transform of the gamma distribution.
s = sw.*(1 + xps.b.*mu2./d_iso).^(-1*(d_iso.^2./mu2));        

