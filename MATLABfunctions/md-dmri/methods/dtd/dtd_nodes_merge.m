function dtd_nodes = dtd_nodes_merge(dtd_nodes1,dtd_nodes2)

[n1,par1,perp1,theta1,phi1] = dtd_nodes2par(dtd_nodes1);
[n2,par2,perp2,theta2,phi2] = dtd_nodes2par(dtd_nodes2);

n = n1 + n2;
par = [par1; par2];
perp = [perp1; perp2];
theta = [theta1; theta2];
phi = [phi1; phi2];

dtd_nodes = [par'; perp'; theta'; phi'];
dtd_nodes = [n; dtd_nodes(:)];
