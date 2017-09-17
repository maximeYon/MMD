function q = gwf_q_from_g(g, dt)
% function q = gwf_q_from_g(g, dt)
%
% Compute q(t) from gradient g(t), assuming gamma for the 1H nuclei
%
% Input
% g  - gradient trajectory
% dt - time step
%
% Output
% q  - q-trajectory

q = cumsum(msf_const_gamma('1H') * g * dt);


