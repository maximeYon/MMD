function dtor1r2d_nodes = dtor1r2d_rand(n,dmin,dmax,rmin,rmax,r1min,r1max,r2min,r2max,opt)

par = dmin*(dmax/dmin).^rand(1,n);
perp = dmin*(dmax/dmin).^rand(1,n);
theta = acos(2*rand(1,n)-1);
phi = 2*pi*rand(1,n);
d0 = dmin*(dmax/dmin).^rand(1,n);
rpar = rmin*(rmax/rmin).^rand(1,n);
rperp = rmin*(rmax/rmin).^rand(1,n);
r1 = r1min*(r1max/r1min).^rand(1,n);
r2 = r2min*(r2max/r2min).^rand(1,n);

if opt.dtor1r2d.no_ivim
    ind = d0<par;
    d0(ind) = par(ind);
    ind = d0<perp;
    d0(ind) = perp(ind);
end

% [par' perp' theta' phi' d0' rpar' rperp']

dtor1r2d_nodes = [par; perp; theta; phi; d0; rpar; rperp; r1; r2];
dtor1r2d_nodes = [n; dtor1r2d_nodes(:)];


