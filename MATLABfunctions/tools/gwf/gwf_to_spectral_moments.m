function [M0,M1] = gwf_to_spectral_moments(gwf, rf, dt)
% function [M0,M1] = gwf_to_spectral_moments(gwf, rf, dt)
%

if (isempty(rf)), rf = gwf_to_rf(gwf); end

g_eff = gwf .* repmat(rf, 1, 3);
q_eff = cumsum(g_eff, 1) * dt;

% M1 referred to as bV_omega in Nilsson et al (2017) NMR Biomed, Eq. 24
M0 = msf_const_gamma('1H')^2 * (q_eff') * q_eff * dt;
M1 = msf_const_gamma('1H')^2 * (g_eff') * g_eff * dt;
