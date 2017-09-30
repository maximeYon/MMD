function q = gwf_to_q(gwf, rf, dt)
% function q = gwf_to_q(gwf, rf, dt)

if (isempty(rf)), rf = gwf_to_rf(gwf); end

g_eff = gwf .* repmat(rf, 1, 3);

q = msf_const_gamma('1H') * cumsum( g_eff, 1 ) * dt;
