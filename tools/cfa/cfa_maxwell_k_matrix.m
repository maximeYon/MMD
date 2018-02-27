function k = cfa_maxwell_k_matrix(gwf, rf, dt, B0)
% function k = cfa_maxwell_k_matrix(gwf, rf, dt, B0)
% 
% k    - matrix that determines size of the effect of maxwell terms.
% gwf  - gradient waveform, size n x 3 
% rf   - effect of gradient waveform, n x 1 
% 
% References:
% Bernstein et al., Concomitant gradient terms in phase contrast MR: 
% analysis and correction. Magn Reson Med, 1998. 39(2): p. 300-8.
% 
% Baron et al., The effect of concomitant gradient fields on diffusion 
% tensor imaging. Magn Reson Med, 2012. 68(4): p. 1190-201.
%
% Szczepankiewicz and Nilsson, Maxwell-compensated waveform design for 
% asymmetric diffusion encoding. ISMRM, 2018, Paris, France
% Download abstract at: https://goo.gl/vVGQq2

Gx = gwf(:,1);
Gy = gwf(:,2);
Gz = gwf(:,3);

t = [
       sum(Gz.*Gz.*rf),                     0,                      -2*sum(Gx.*Gz.*rf)  ;
                     0,       sum(Gz.*Gz.*rf),                      -2*sum(Gy.*Gz.*rf)  ;
    -2*sum(Gx.*Gz.*rf),    -2*sum(Gy.*Gz.*rf),    4*(sum(Gx.*Gx.*rf) + sum(Gy.*Gy.*rf)) ;
    ];

k = msf_const_gamma / 2 / pi * t / (4 * B0) * dt;

