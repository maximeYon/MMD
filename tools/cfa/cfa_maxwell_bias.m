function [c, c_s, c_p] = cfa_maxwell_bias(gwf, rf, dt, ips)
% function [c, c_s, c_p] = cfa_maxwell_bias(gwf, rf, dt, ips)
%
% gwf - gradient waveform of size n x 3
% rf  - effect of rf, size 1 x 3
% dt  - time step
% ips - imaging parameter structure (see function "csf_ips_example")
%
% Baron et al., The effect of concomitant gradient fields on diffusion 
% tensor imaging. Magn Reson Med, 2012. 68(4): p. 1190-201.

k   = cfa_maxwell_k_matrix(gwf, rf, dt, ips.B0);

c_s = cfa_maxwell_bias_slice(k, ips.o(3,:), ips.r_xyz, ips.res(3));
c_p = cfa_maxwell_bias_phase(k, ips.o(2,:), ips.r_xyz, ips.T2s, ips.kpv);

% Total bias is the product of phase and slice biases
c   = c_s .* c_p;



