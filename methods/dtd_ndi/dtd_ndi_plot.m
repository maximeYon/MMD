function dtd_ndi_plot(S, xps, h, h2)
% function dtd_ndi_plot(S, xps, h, h2)

if (nargin < 4), h2 = []; end

ind = xps.b_delta > 0.95;

m = mplot_s_vs_b_by_b_delta(S, xps, ...
    @dtd_ndi_1d_data2fit, @dtd_ndi_1d_fit2data, @dtd_ndi_opt, h, h2, ind);

title(h, sprintf('v(neurite) = %0.2f, v(csf) = %0.2f', m(2), m(3)));