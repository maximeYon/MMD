function dtor1r2d = dtor1r2d_darwin(stemp, xps, dtor1r2d, opt)

n_nodes = opt.dtor1r2d.n_in;
n_darwin = opt.dtor1r2d.n_darwin;

for niter = 1:n_darwin    
    n_in = dtor1r2d(1);
    n_max = min([(n_in-opt.dtor1r2d.n_kill) opt.dtor1r2d.n_out dtor1r2d(1)]);

    dtor1r2d_nodes1 = dtor1r2d_dist2nodes(dtor1r2d);
    ind = 1:n_max;
    dtor1r2d_nodes1 = dtor1r2d_nodes_select(dtor1r2d_nodes1,ind);
    ind = 1 + floor((n_max-1)*linspace(0,1,n_nodes).^3);
    ind(ind<1) = 1;
    %dtor1r2d_nodes1(1)
    %ind
    if dtor1r2d_nodes1(1) == 0
        dtor1r2d = [];
        break
    end
    
    dtor1r2d_nodes2 = dtor1r2d_nodes_select(dtor1r2d_nodes1,ind);
    dtor1r2d_nodes2 = dtor1r2d_nodes_mutate(dtor1r2d_nodes2,opt);

    dtor1r2d_nodes = dtor1r2d_nodes_merge(dtor1r2d_nodes1,dtor1r2d_nodes2);

    dtor1r2d = dtor1r2d_data2dtor1r2d(stemp, xps, dtor1r2d_nodes);
end

if ~isempty(dtor1r2d)
    dtor1r2d_nodes = dtor1r2d_dist2nodes(dtor1r2d);
    n_max = min([opt.dtor1r2d.n_out dtor1r2d(1)]);
    ind = 1:n_max;
    dtor1r2d_nodes = dtor1r2d_nodes_select(dtor1r2d_nodes,ind);
    dtor1r2d = dtor1r2d_data2dtor1r2d(stemp, xps, dtor1r2d_nodes);
end
