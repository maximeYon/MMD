function dtod_nodes = dtod_nodes_merge(dtod_nodes1,dtod_nodes2)

[n1,par1,perp1,theta1,phi1,d01,rpar1,rperp1] = dtod_nodes2par(dtod_nodes1);
[n2,par2,perp2,theta2,phi2,d02,rpar2,rperp2] = dtod_nodes2par(dtod_nodes2);

n = n1 + n2;
par = [par1; par2];
perp = [perp1; perp2];
theta = [theta1; theta2];
phi = [phi1; phi2];
d0 = [d01; d02];
rpar = [rpar1; rpar2];
rperp = [rperp1; rperp2];

dtod_nodes = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'];
dtod_nodes = [n; dtod_nodes(:)];
