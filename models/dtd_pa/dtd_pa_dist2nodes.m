function dtd_nodes = dtd_pa_dist2nodes(dtd)

[n,par,perp,w] = dtd_pa_dist2par(dtd);

dtd_nodes = [par'; perp'];
dtd_nodes = [n; dtd_nodes(:)];
