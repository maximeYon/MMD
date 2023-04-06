function b = gwf_b_from_q(q, dt)
% function b = gwf_b_from_q(q, dt)
%
% Compute the b-value from a q-trajectory
%
% q  - the q-trajectory of size (1, 2, or 3  x  n )
% dt - the time step

b = sum( q.^2 * dt, 2 );


