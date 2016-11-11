function erf_plot(S, xps, h, h2)
% function erf_plot(S, xps, h, h2)

if (nargin < 4), h2 = []; end


ind = (S > 0.1 * max(S));

m = mplot_s_vs_b_by_b_delta(S, xps, ...
    @erf_1d_data2fit, @erf_1d_fit2data, @erf_opt, h, h2, ind);

title(h, sprintf('D_I = %0.2f, D(delta) = %0.2f', m(2) * 1e9, m(3)));
