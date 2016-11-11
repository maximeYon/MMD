function dtd = dtd_pa_nodesw2dist(dtd_nodes,w)

[n,par,perp] = dtd_pa_nodes2par(dtd_nodes);

dtd = [par'; perp'; w'];
dtd = [n; dtd(:)];
