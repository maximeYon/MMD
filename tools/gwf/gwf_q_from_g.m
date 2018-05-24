function q = gwf_q_from_g(gwf, rf, dt)
% function q = gwf_q_from_g(gwf, dt)
%
% Compute q(t) from gradient g(t), assuming gamma for the 1H nuclei
%
% Input
% g  - gradient trajectory as executed on amplifier
% rf - effect of rf waveform
% dt - time step
%
% Output
% q  - q-trajectory

g_eff = gwf_to_g_eff(gwf, rf, dt);

q = cumsum(msf_const_gamma('1H') * g_eff * dt);


