function dtd = dtd_nodesw2dist(dtd_nodes,w)

[n,par,perp,theta,phi] = dtd_nodes2par(dtd_nodes);

dtd = [par'; perp'; theta'; phi'; w'];
dtd = [n; dtd(:)];
