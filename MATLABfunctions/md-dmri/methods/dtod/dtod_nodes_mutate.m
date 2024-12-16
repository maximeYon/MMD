function dtod_nodes_out = dtod_mutate(dtod_nodes,opt)

[n,par,perp,theta,phi,d0,rpar,rperp] = dtod_nodes2par(dtod_nodes);

par = par.*(1+opt.dtod.dfuzz*randn(n,1));
perp = perp.*(1+opt.dtod.dfuzz*randn(n,1));
theta = theta + opt.dtod.ofuzz*randn(n,1);
phi = phi + opt.dtod.ofuzz*randn(n,1);
d0 = d0.*(1+opt.dtod.dfuzz*randn(n,1));
rpar = rpar.*(1+opt.dtod.rfuzz*randn(n,1));
rperp = rperp.*(1+opt.dtod.rfuzz*randn(n,1));

par(par>opt.dtod.dmax) = opt.dtod.dmax;
perp(perp>opt.dtod.dmax) = opt.dtod.dmax;
d0(d0>opt.dtod.dmax) = opt.dtod.dmax;
rpar(rpar>opt.dtod.rmax) = opt.dtod.rmax;
rperp(rperp>opt.dtod.rmax) = opt.dtod.rmax;
par(par<opt.dtod.dmin) = opt.dtod.dmin;
perp(perp<opt.dtod.dmin) = opt.dtod.dmin;
d0(d0<opt.dtod.dmin) = opt.dtod.dmin;
rpar(rpar<opt.dtod.rmin) = opt.dtod.rmin;
rperp(rperp<opt.dtod.rmin) = opt.dtod.rmin;

if opt.dtod.no_ivim
    ind = d0<par;
    d0(ind) = par(ind);
    ind = d0<perp;
    d0(ind) = perp(ind);
end

dtod_nodes_out = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'];
dtod_nodes_out = [n; dtod_nodes_out(:)];
