function dtod_nodes = dtod_rand(n,dmin,dmax,rmin,rmax)

par = dmin*(dmax/dmin).^rand(1,n);
perp = dmin*(dmax/dmin).^rand(1,n);
theta = acos(2*rand(1,n)-1);
phi = 2*pi*rand(1,n);
d0 = dmin*(dmax/dmin).^rand(1,n);
rpar = rmin*(rmax/rmin).^rand(1,n);
rperp = rmin*(rmax/rmin).^rand(1,n);

ind = d0<par;
d0(ind) = par(ind);
ind = d0<perp;
d0(ind) = perp(ind);

% [par' perp' theta' phi' d0' rpar' rperp']

dtod_nodes = [par; perp; theta; phi; d0; rpar; rperp];
dtod_nodes = [n; dtod_nodes(:)];


