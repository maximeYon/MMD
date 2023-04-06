function dtod = dtod_nodesw2dist(dtod_nodes,w)

[n,par,perp,theta,phi,d0,rpar,rperp] = dtod_nodes2par(dtod_nodes);

dtod = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; w'];
dtod = [n; dtod(:)];
