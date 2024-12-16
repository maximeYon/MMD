function dtod = dtod_proliferation(s, xps, opt)

dmin = opt.dtod.dmin;
dmax = opt.dtod.dmax;
rmin = opt.dtod.rmin;
rmax = opt.dtod.rmax;
n_nodes = opt.dtod.n_in;
n_proliferation = opt.dtod.n_proliferation;

dtod_nodes1 = [];
for niter = 1:n_proliferation    
    dtod_nodes2 = dtod_rand(n_nodes,dmin,dmax,rmin,rmax,opt);
    dtod_nodes = dtod_nodes_merge(dtod_nodes1,dtod_nodes2);

    dtod = dtod_data2dtod(s, xps, dtod_nodes);
    dtod_nodes1 = dtod_dist2nodes(dtod);

%     m = dtod_dtod2m(dtod,opt);
%     [n,dpar,dperp,theta,phi,d0,rpar,rperp,w] = dtod_dist2par(dtod);
%     [dpar dperp theta phi d0 rpar rperp w]
%     s_fit = dtod_1d_fit2data(m, xps);
%     figure(1), clf
%     plot(1:xps.n,s,'o',1:xps.n,s_fit,'.')
%     title(num2str(niter))
%     pause(1)

end

