function bt = gwf_to_bt(gwf, rf, dt)
% function bt = gwf_to_bt(gwf, rf, dt)
%
% gwf - gradient waveform of size N x 3
% rf  - effect of rf pulses (range -1 to 1), size N x 1
% dt  - time step of waveform
%
% Following notation in Westin et al (2016) NeuroImage 135

q = gwf_to_q(gwf, rf, dt);

bt = q' * q * dt;
