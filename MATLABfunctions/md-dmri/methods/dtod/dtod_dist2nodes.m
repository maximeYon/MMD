function dtod_nodes = dtod_dist2nodes(dtod)

[n,par,perp,theta,phi,d0,rpar,rperp,w] = dtod_dist2par(dtod);

dtod_nodes = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'];
dtod_nodes = [n; dtod_nodes(:)];
