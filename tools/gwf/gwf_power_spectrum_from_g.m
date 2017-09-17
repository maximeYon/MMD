function [ps,f,df] = gwf_power_spectrum_from_g(g, dt)
% function [ps,f,df] = gwf_power_spectrum_from_g(g, dt)
%
% Compute the power spectrum from a gradient trajectory
%
% Input
% g  - the gradient trajectory
% dt - the time step
%
% Output
% ps - encoding power spectrum
% f  - frequency axis
% df - frequency step


if (size(g,1) ~= 1), error('gradient should be of size 1 x n'); end


[ps,f,df] = gwf_power_spectrum_from_q(gwf_q_from_g(g, dt),dt);
