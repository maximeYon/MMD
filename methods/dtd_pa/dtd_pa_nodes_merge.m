function dtd_nodes = dtd_pa_nodes_merge(dtd_nodes1,dtd_nodes2)

[n1,par1,perp1] = dtd_pa_nodes2par(dtd_nodes1);
[n2,par2,perp2] = dtd_pa_nodes2par(dtd_nodes2);

n = n1 + n2;
par = [par1; par2];
perp = [perp1; perp2];

dtd_nodes = [par'; perp'];
dtd_nodes = [n; dtd_nodes(:)];
