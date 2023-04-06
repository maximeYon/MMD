function s = dti_nls_1d_fit2data(m, xps)
% function s = dti_nls_1d_fit2data(m, xps)
%
% Predict the signal 's' as a 1 x xps.n vector from the model parameters
% in 'm'. 

dt = zeros(1,6);
dt(1:6) = m(2:7);

s = m(1) * exp(-dt * xps.bt')';
%                1x6  (n x 6)' = (1x6 * 6 x n) = 1 x n