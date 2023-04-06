function dtd_nodes_out = dtd_pa_nodes_select(dtd_nodes,ind)

[n,par,perp] = dtd_pa_nodes2par(dtd_nodes);

n = numel(ind);
par = par(ind);
perp = perp(ind);

dtd_nodes_out = [par'; perp'];
dtd_nodes_out = [n; dtd_nodes_out(:)];
