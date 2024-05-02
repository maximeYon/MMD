function dtd_nodes = dtd_dist2nodes(dtd)

[n,par,perp,theta,phi,w] = dtd_dist2par(dtd);

dtd_nodes = [par'; perp'; theta'; phi'];
dtd_nodes = [n; dtd_nodes(:)];
