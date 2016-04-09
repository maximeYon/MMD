function dtd_nodes_out = dtd_mutate(dtd_nodes,opt)

[n,par,perp,theta,phi] = dtd_nodes2par(dtd_nodes);

par = par.*(1+opt.dtd.dfuzz*randn(n,1));
perp = perp.*(1+opt.dtd.dfuzz*randn(n,1));
theta = theta + opt.dtd.ofuzz*randn(n,1);
phi = phi + opt.dtd.ofuzz*randn(n,1);

par(par>opt.dtd.dmax) = opt.dtd.dmax;
perp(perp>opt.dtd.dmax) = opt.dtd.dmax;
par(par<opt.dtd.dmin) = opt.dtd.dmin;
perp(perp<opt.dtd.dmin) = opt.dtd.dmin;

dtd_nodes_out = [par'; perp'; theta'; phi'];
dtd_nodes_out = [n; dtd_nodes_out(:)];
