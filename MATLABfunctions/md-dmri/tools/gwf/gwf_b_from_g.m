function b = gwf_b_from_g(gwf,rf,dt)
% function b = gwf_b_from_g(gwf,rf,dt)
%
% Compute the b-value from a gradient waveform
%
% g  - gradient waveform as executed on amplifier
% rf - effect of rf waveform
% dt - time step

b = gwf_b_from_q(gwf_q_from_g(gwf, rf, dt), dt);


