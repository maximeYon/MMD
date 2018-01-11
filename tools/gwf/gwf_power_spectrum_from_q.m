function [ps,f,df] = gwf_power_spectrum_from_q(q, dt)
% function [ps,f,df] = gwf_power_spectrum_from_q(q, dt)
%
% Compute the encoding power spectrum from q(t) 
%
% Input
% q  - the q-trajectory
% dt - time step
%
% Output
% ps - encoding power spectrum
% f  - frequency axis
% df - frequency step

ps = abs(fftshift(fft(q * dt))).^2;
f  = linspace(-1/dt, 1/dt, numel(ps)) / 2;
df = f(2) - f(1);

