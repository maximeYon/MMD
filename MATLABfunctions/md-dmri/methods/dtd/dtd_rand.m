function dtd_nodes = dtd_rand(n,dmin,dmax)

par = dmin*(dmax/dmin).^rand(1,n);
perp = dmin*(dmax/dmin).^rand(1,n);
theta = acos(2*rand(1,n)-1);
phi = 2*pi*rand(1,n);

dtd_nodes = [par; perp; theta; phi];
dtd_nodes = [n; dtd_nodes(:)];


