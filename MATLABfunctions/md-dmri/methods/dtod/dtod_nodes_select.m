function dtod_nodes_out = dtod_nodes_select(dtod_nodes,ind)

[n,par,perp,theta,phi,d0,rpar,rperp] = dtod_nodes2par(dtod_nodes);

n = numel(ind);
par = par(ind);
perp = perp(ind);
theta = theta(ind);
phi = phi(ind);
d0 = d0(ind);
rpar = rpar(ind);
rperp = rperp(ind);

dtod_nodes_out = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'];
dtod_nodes_out = [n; dtod_nodes_out(:)];
