function dtd = dtd_extinction(stemp, bt_mx6, dtd, opt)

n_nodes = opt.dtd.n_in;
n_extinction = opt.dtd.n_extinction;

for niter = 1:n_extinction    
    n_in = dtd(1);
    n_max = min([(n_in-opt.dtd.n_kill) opt.dtd.n_out dtd(1)]);

    dtd_nodes1 = dtd_dist2nodes(dtd);
    ind = 1:n_max;
    dtd_nodes1 = dtd_nodes_select(dtd_nodes1,ind);
    ind = 1 + floor((n_max-1)*linspace(0,1,n_nodes).^3);
    ind(ind<1) = 1;
    %dtd_nodes1(1)
    %ind
    if dtd_nodes1(1) == 0
        dtd = [];
        break
    end
    
    dtd_nodes2 = dtd_nodes_select(dtd_nodes1,ind);
    dtd_nodes2 = dtd_nodes_mutate(dtd_nodes2,opt);

    dtd_nodes = dtd_nodes_merge(dtd_nodes1,dtd_nodes2);

    dtd = dtd_data2dtd(stemp,bt_mx6,dtd_nodes);    
end

if ~isempty(dtd)
    dtd_nodes = dtd_dist2nodes(dtd);
    n_max = min([opt.dtd.n_out dtd(1)]);
    ind = 1:n_max;
    dtd_nodes = dtd_nodes_select(dtd_nodes,ind);
    dtd = dtd_data2dtd(stemp,bt_mx6,dtd_nodes);    
end
