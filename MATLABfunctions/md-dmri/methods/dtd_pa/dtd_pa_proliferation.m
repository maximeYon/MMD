function dtd = dtd_pa_proliferation(stemp, b, b_delta, opt)

dmin = opt.dtd_pa.dmin;
dmax = opt.dtd_pa.dmax;
n_nodes = opt.dtd_pa.n_in;
n_proliferation = opt.dtd_pa.n_proliferation;

dtd_nodes1 = [];
for niter = 1:n_proliferation    
    dtd_nodes2 = dtd_pa_rand(n_nodes,dmin,dmax);
    dtd_nodes = dtd_pa_nodes_merge(dtd_nodes1,dtd_nodes2);

    dtd = dtd_pa_data2dtd(stemp, b, b_delta, dtd_nodes);    
    dtd_nodes1 = dtd_pa_dist2nodes(dtd);
end

