function dtd_gamma_plot(S, xps, h, h2)
% function dtd_gamma_plot(S, xps, h, h2)

if (nargin < 4), h2 = []; end

if any(~isreal(S(:)))
    S = abs(S);
end

m = mplot_s_vs_b_by_b_delta(S, xps, ...
    @dtd_gamma_1d_data2fit, @dtd_gamma_1d_fit2data, @dtd_gamma_opt, h, h2);

s0 = m(1);
md = m(2) / 1e-9;
mu2iso_n = m(3)/(m(2)^2);
mu2aniso_n = m(4)/(m(2)^2);

Vl = 5/2 * m(4);

uFA = sqrt(3/2) * sqrt( Vl ./ (Vl + m(3) + m(2).^2) );

title_str = {...
    ['MD = ' num2str(md, 2) ' ',char(181),'m^2/ms']; ...
    ['\mu_2^{iso}/MD^2 = ' num2str(mu2iso_n, 2) '   \mu_2^{aniso}/MD^2 = ' num2str(mu2aniso_n, 2)]};
    %['\mu_2^{iso}/MD^2 = ' num2str(mu2iso_n, 2) ' \mu_2^{aniso}/MD^2 = ' num2str(mu2aniso_n, 2) ' ',char(181),'FA = ' num2str(uFA, 2)]};

title(h, title_str)