function gamma_plot(S, xps, h, h2)
% function gamma_plot(S, xps, h, h2)

if (nargin < 4), h2 = []; end

mplot_s_vs_b_by_b_delta(S, xps, ...
    @gamma_1d_data2fit, @gamma_1d_fit2data, @gamma_opt, h, h2);