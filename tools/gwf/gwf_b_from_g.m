function b = gwf_b_from_g(g,dt)
% function b = gwf_b_from_g(g,dt)
%
% Compute the b-value from a gradient waveform
%
% g  - gradient waveform in (1, 2 or 3   x  n )
% dt - time step


b = gwf_b_from_q(gwf_q_from_g(g, dt), dt);


