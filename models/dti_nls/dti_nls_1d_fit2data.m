function s = dti_nls_1d_fit2data(m, xps)
% function s = dti_nls_1d_fit2data(m, xps)

C = [...
    m(2) m(5) m(7);
    0    m(3) m(6);
    0     0   m(4)];

dt = dtd_3x3_to_1x6(C' * C);

s = m(1) * exp(-(dt * xps.bt'))';
%                1x6  (n x 6)' = (1x6 * 6 x n) = 1 x n