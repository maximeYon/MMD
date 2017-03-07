function s = dki_pa_1d_fit2data(m, xps)
% function s = dki_pa_1d_fit2data(m, xps)
%
% m(1) - s0
% m(2) - md
% m(3) - v(d)

% define data from fitted parameters
s = m(1) * exp(-xps.b * m(2) + 0.5 * xps.b.^2 * m(3));
