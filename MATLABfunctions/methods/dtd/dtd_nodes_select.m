function dtd_nodes_out = dtd_nodes_select(dtd_nodes,ind)

[n,par,perp,theta,phi] = dtd_nodes2par(dtd_nodes);

n = numel(ind);
par = par(ind);
perp = perp(ind);
theta = theta(ind);
phi = phi(ind);

dtd_nodes_out = [par'; perp'; theta'; phi'];
dtd_nodes_out = [n; dtd_nodes_out(:)];
