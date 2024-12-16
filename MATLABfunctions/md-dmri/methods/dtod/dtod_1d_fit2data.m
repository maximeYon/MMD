function s = dtod_1d_fit2data(m, xps)

if m(1)>0
    dtod = dtod_m2dtod(m);
    [dtod_nx6,w] = dtod_dist2nx6w(dtod,xps);

%     figure(3), clf
%     subplot(3,1,1)
%     plot(1:size(xps.bto,2),xps.bto(2637,:),'-')
%     subplot(3,1,2)
%     plot(1:size(xps.bto,2),dtod_nx6,'-')
    %[dtod_nx6,w] = dtod_dist2nx6w(dtod);
    k = exp(-xps.btomega*dtod_nx6');
%     k = exp(-xps.bt*dtod_nx6').*(1-exp(-xps.tr*r1')).*exp(-xps.te*r2');
    s = k*w;
else
    s = zeros(xps.n,1);
end

