function m = dtd_pa_1d_data2fit(signal, xps, opt, ind)
% function m = dtd_pa_1d_data2fit(signal, xps, opt, ind)
%
% Size-shape diffusion tensor distribution
% Assumes powder averaging
% de Almeida Martins and Topgaard, Phys. Rev. Lett. 116, 087601 (2016).
% http://dx.doi.org/10.1103/PhysRevLett.116.087601
% Modified for general b-tensor shapes

if (nargin < 4), ind = ones(size(signal)) > 0; end

b = xps.b(ind,:);
b_delta = xps.b_delta(ind,:);
stemp = signal(ind);

dtd = dtd_pa_proliferation(stemp, b, b_delta, opt);

dtd = dtd_pa_extinction(stemp, b, b_delta, dtd, opt);

m = dtd_pa_dtd2m(dtd,opt);

if (opt.dtd_pa.do_plot)
    figure(1), clf
    signal_fit = dtd_pa_1d_fit2data(m, xps);
    plot(1:xps.n,signal,'o',1:xps.n,signal_fit,'x');
    pause(0.05);
end