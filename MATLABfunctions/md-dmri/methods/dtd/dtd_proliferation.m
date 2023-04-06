function dtd = dtd_proliferation(stemp, bt_mx6, opt)

dmin = opt.dtd.dmin;
dmax = opt.dtd.dmax;
n_nodes = opt.dtd.n_in;
n_proliferation = opt.dtd.n_proliferation;

dtd_nodes1 = [];
for niter = 1:n_proliferation    
    dtd_nodes2 = dtd_rand(n_nodes,dmin,dmax);
    dtd_nodes = dtd_nodes_merge(dtd_nodes1,dtd_nodes2);

    dtd = dtd_data2dtd(stemp,bt_mx6,dtd_nodes);    
    dtd_nodes1 = dtd_dist2nodes(dtd);
end

