function dtod = dtod_darwin(stemp, xps, dtod, opt)

n_nodes = opt.dtod.n_in;
n_darwin = opt.dtod.n_darwin;

for niter = 1:n_darwin    
    n_in = dtod(1);
    n_max = min([(n_in-opt.dtod.n_kill) opt.dtod.n_out dtod(1)]);

    dtod_nodes1 = dtod_dist2nodes(dtod);
    ind = 1:n_max;
    dtod_nodes1 = dtod_nodes_select(dtod_nodes1,ind);
    ind = 1 + floor((n_max-1)*linspace(0,1,n_nodes).^3);
    ind(ind<1) = 1;
    %dtod_nodes1(1)
    %ind
    if dtod_nodes1(1) == 0
        dtod = [];
        break
    end
    
    dtod_nodes2 = dtod_nodes_select(dtod_nodes1,ind);
    dtod_nodes2 = dtod_nodes_mutate(dtod_nodes2,opt);

    dtod_nodes = dtod_nodes_merge(dtod_nodes1,dtod_nodes2);

    dtod = dtod_data2dtod(stemp, xps, dtod_nodes);
end

if ~isempty(dtod)
    dtod_nodes = dtod_dist2nodes(dtod);
    n_max = min([opt.dtod.n_out dtod(1)]);
    ind = 1:n_max;
    dtod_nodes = dtod_nodes_select(dtod_nodes,ind);
    dtod = dtod_data2dtod(stemp, xps, dtod_nodes);
end
