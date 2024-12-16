function dtor1r2d = dtor1r2d_proliferation(s, xps, opt)

dmin = opt.dtor1r2d.dmin;
dmax = opt.dtor1r2d.dmax;
rmin = opt.dtor1r2d.rmin;
rmax = opt.dtor1r2d.rmax;
r1min = opt.dtor1r2d.r1min;
r1max = opt.dtor1r2d.r1max;
r2min = opt.dtor1r2d.r2min;
r2max = opt.dtor1r2d.r2max;
n_nodes = opt.dtor1r2d.n_in;
n_proliferation = opt.dtor1r2d.n_proliferation;

dtor1r2d_nodes1 = [];
for niter = 1:n_proliferation    
    dtor1r2d_nodes2 = dtor1r2d_rand(n_nodes,dmin,dmax,rmin,rmax,r1min,r1max,r2min,r2max,opt);
    dtor1r2d_nodes = dtor1r2d_nodes_merge(dtor1r2d_nodes1,dtor1r2d_nodes2);

    dtor1r2d = dtor1r2d_data2dtor1r2d(s, xps, dtor1r2d_nodes);
    dtor1r2d_nodes1 = dtor1r2d_dist2nodes(dtor1r2d);

%     m = dtor1r2d_dtor1r2d2m(dtor1r2d,opt);
%     [n,dpar,dperp,theta,phi,d0,rpar,rperp,w] = dtor1r2d_dist2par(dtor1r2d);
%     [dpar dperp theta phi d0 rpar rperp w]
%     s_fit = dtor1r2d_1d_fit2data(m, xps);
%     figure(1), clf
%     plot(1:xps.n,s,'o',1:xps.n,s_fit,'.')
%     title(num2str(niter))
%     pause(1)

end

