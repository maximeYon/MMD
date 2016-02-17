function m = dti_nls_1d_data2fit(signal, xps, opt)
% function m = dti_nls_1d_data2fit(signal, xps, opt)
%
% Yields a 1x7 vector 'm' with the fit parameters
% m(1)   - s0
% m(2:7) - cholesky decomposition of the diffusion tensor

signal = double(signal);

unit_to_SI = [max(signal) [1 1 1 1 1 1] * 10^-(9/2)];

fun = @(m,e) dti_nls_1d_fit2data(m .* unit_to_SI, xps);

m_guess = [1 1 1 1  0  0  0];
m_lb    = [0 0 0 0 -9 -9 -9];
m_ub    = [2 9 9 9 +9 +9 +9];

m = lsqcurvefit(fun, m_guess, [], signal, ...
    m_lb, m_ub, opt.dti_nls.lsq_opts);

m = m .* unit_to_SI;