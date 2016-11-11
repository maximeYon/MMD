function s = fexi11_1d_fit2data(m, xps)
% function s = fexi11_1d_fit2data(m, xps)
%

ind_non_eq = xps.mde_b1 > 0; % index of filtered data

ADC0  = m(1);
sigma = m(2);
AXR   = m(3);
S0    = m(4:end);

ADC = ADC0 * (1 - ind_non_eq .* sigma .* exp(-AXR * xps.mde_tm12));
s   = exp(-ADC .* xps.mde_b2);

s   = S0(xps.s_ind)' .* s;

