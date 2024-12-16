function dtor1r2d_nodes = dtor1r2d_nodes_merge(dtor1r2d_nodes1,dtor1r2d_nodes2)

[n1,par1,perp1,theta1,phi1,d01,rpar1,rperp1,r11,r21] = dtor1r2d_nodes2par(dtor1r2d_nodes1);
[n2,par2,perp2,theta2,phi2,d02,rpar2,rperp2,r12,r22] = dtor1r2d_nodes2par(dtor1r2d_nodes2);

n = n1 + n2;
par = [par1; par2];
perp = [perp1; perp2];
theta = [theta1; theta2];
phi = [phi1; phi2];
d0 = [d01; d02];
rpar = [rpar1; rpar2];
rperp = [rperp1; rperp2];
r1 = [r11; r12];
r2 = [r21; r22];

dtor1r2d_nodes = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; r1'; r2'];
dtor1r2d_nodes = [n; dtor1r2d_nodes(:)];
