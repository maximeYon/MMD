function erf_plot(S, xps, h, h2)
% function erf_plot(S, xps, h, h2)

if (nargin < 4), h2 = []; end

mplot_s_vs_b_by_b_delta(S, xps, ...
    @erf_1d_data2fit, @erf_1d_fit2data, @erf_opt, h, h2);