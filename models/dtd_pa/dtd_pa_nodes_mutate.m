function dtd_nodes_out = dtd_pa_mutate(dtd_nodes,opt)

[n,par,perp] = dtd_pa_nodes2par(dtd_nodes);

par = par.*(1+opt.dtd_pa.dfuzz*randn(n,1));
perp = perp.*(1+opt.dtd_pa.dfuzz*randn(n,1));

par(par>opt.dtd_pa.dmax) = opt.dtd_pa.dmax;
perp(perp>opt.dtd_pa.dmax) = opt.dtd_pa.dmax;
par(par<opt.dtd_pa.dmin) = opt.dtd_pa.dmin;
perp(perp<opt.dtd_pa.dmin) = opt.dtd_pa.dmin;

dtd_nodes_out = [par'; perp'];
dtd_nodes_out = [n; dtd_nodes_out(:)];
