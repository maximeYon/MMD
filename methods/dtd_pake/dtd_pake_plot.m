function dtd_pake_plot(S, xps, h, h2)
% function dtd_pake_plot(S, xps, h, h2)

if (nargin < 4), h2 = []; end

if any(~isreal(S(:)))
    S = abs(S);
end

m = mplot_s_vs_b_by_b_delta(S, xps, ...
    @dtd_pake_1d_data2fit, @dtd_pake_1d_fit2data, @dtd_pake_opt, h, h2);

s0 = m(1);
d_iso = m(2);
d_delta = m(3);

title_str = {...
    ['MD = ' num2str(d_iso/1e-9, 2) ' ',char(181),'m^2/ms']; ...
    ['\itD\rm_{\Delta} = ' num2str(d_delta, 2) ]};

title(h, title_str)