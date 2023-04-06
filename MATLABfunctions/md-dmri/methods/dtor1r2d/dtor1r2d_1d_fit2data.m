function s = dtor1r2d_1d_fit2data(m, xps)

if m(1)>0
    dtor1r2d = dtor1r2d_m2dtor1r2d(m);
    [dtor1r2d_nx6,r1,r2,w] = dtor1r2d_dist2nx6r1r2w(dtor1r2d,xps);

%     figure(3), clf
%     subplot(3,1,1)
%     plot(1:size(xps.bto,2),xps.bto(2637,:),'-')
%     subplot(3,1,2)
%     plot(1:size(xps.bto,2),dtor1r2d_nx6,'-')
    %[dtor1r2d_nx6,w] = dtor1r2d_dist2nx6w(dtor1r2d);
%     k = exp(-xps.btomega*dtor1r2d_nx6');
    k = exp(-xps.btomega*dtor1r2d_nx6').*(1-exp(-xps.tr*r1')).*exp(-xps.te*r2');
%     k = exp(-xps.bt*dtor1r2d_nx6').*(1-exp(-xps.tr*r1')).*exp(-xps.te*r2');
    s = k*w;
else
    s = zeros(xps.n,1);
end

