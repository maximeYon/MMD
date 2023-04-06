function s = dtd_1d_fit2data(m, xps)

if m(1)>0
    dtd = dtd_m2dtd(m);
    [dtd_nx6,w] = dtd_dist2nx6w(dtd);
    k = exp(-xps.bt*dtd_nx6');
    s = k*w;
else
    s = zeros(xps.n,1);
end

