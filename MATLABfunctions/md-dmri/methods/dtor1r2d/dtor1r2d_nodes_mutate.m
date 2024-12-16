function dtor1r2d_nodes_out = dtor1r2d_mutate(dtor1r2d_nodes,opt)

[n,par,perp,theta,phi,d0,rpar,rperp,r1,r2] = dtor1r2d_nodes2par(dtor1r2d_nodes);

par = par.*(1+opt.dtor1r2d.dfuzz*randn(n,1));
perp = perp.*(1+opt.dtor1r2d.dfuzz*randn(n,1));
theta = theta + opt.dtor1r2d.ofuzz*randn(n,1);
phi = phi + opt.dtor1r2d.ofuzz*randn(n,1);
d0 = d0.*(1+opt.dtor1r2d.dfuzz*randn(n,1));
rpar = rpar.*(1+opt.dtor1r2d.rfuzz*randn(n,1));
rperp = rperp.*(1+opt.dtor1r2d.rfuzz*randn(n,1));
r1 = r1.*(1+opt.dtor1r2d.r1fuzz*randn(n,1));
r2 = r2.*(1+opt.dtor1r2d.r2fuzz*randn(n,1));

par(par>opt.dtor1r2d.dmax) = opt.dtor1r2d.dmax;
perp(perp>opt.dtor1r2d.dmax) = opt.dtor1r2d.dmax;
d0(d0>opt.dtor1r2d.dmax) = opt.dtor1r2d.dmax;
rpar(rpar>opt.dtor1r2d.rmax) = opt.dtor1r2d.rmax;
rperp(rperp>opt.dtor1r2d.rmax) = opt.dtor1r2d.rmax;
r1(r1>opt.dtor1r2d.r1max) = opt.dtor1r2d.r1max;
r2(r2>opt.dtor1r2d.r2max) = opt.dtor1r2d.r2max;
par(par<opt.dtor1r2d.dmin) = opt.dtor1r2d.dmin;
perp(perp<opt.dtor1r2d.dmin) = opt.dtor1r2d.dmin;
d0(d0<opt.dtor1r2d.dmin) = opt.dtor1r2d.dmin;
rpar(rpar<opt.dtor1r2d.rmin) = opt.dtor1r2d.rmin;
rperp(rperp<opt.dtor1r2d.rmin) = opt.dtor1r2d.rmin;
r1(r1<opt.dtor1r2d.r1min) = opt.dtor1r2d.r1min;
r2(r2<opt.dtor1r2d.r2min) = opt.dtor1r2d.r2min;

if opt.dtor1r2d.no_ivim
    ind = d0<par;
    d0(ind) = par(ind);
    ind = d0<perp;
    d0(ind) = perp(ind);
end

dtor1r2d_nodes_out = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; r1'; r2'];
dtor1r2d_nodes_out = [n; dtor1r2d_nodes_out(:)];
