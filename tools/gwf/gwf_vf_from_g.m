function vf = gwf_vf_from_g(g, dt)
% function vf = gwf_vf_from_g(g, dt)
%
% Input
% g  - gradient trajectory
% dt - time step
%
% Output
% vf - variance (second momemnt) of the encoding power spectrum

vf = sum( (msf_const_gamma('1H') * g).^2 * dt ) ...
    / gwf_b_from_g(g,dt) / (2 * pi)^2;
