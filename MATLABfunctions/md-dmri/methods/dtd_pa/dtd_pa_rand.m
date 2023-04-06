function dtd_nodes = dtd_pa_rand(n,dmin,dmax)

par = dmin*(dmax/dmin).^rand(1,n);
perp = dmin*(dmax/dmin).^rand(1,n);

dtd_nodes = [par; perp];
dtd_nodes = [n; dtd_nodes(:)];


