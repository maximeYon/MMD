function dtd_gamma_plot(S, xps, h, h2)
% function dtd_gamma_plot(S, xps, h, h2)

if (nargin < 4), h2 = []; end

mplot_s_vs_b_by_b_delta(S, xps, ...
    @dtd_gamma_1d_data2fit, @dtd_gamma_1d_fit2data, @dtd_gamma_opt, h, h2);