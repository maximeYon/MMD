function dtod = dtod_merge(dtod1,dtod2)

[n1,par1,perp1,theta1,phi1,d01,rpar1,rperp1,w1] = dtod_dist2par(dtod1);
[n2,par2,perp2,theta2,phi2,d02,rpar2,rperp2,w2] = dtod_dist2par(dtod2);

n = n1 + n2;
par = [par1; par2];
perp = [perp1; perp2];
theta = [theta1; theta2];
phi = [phi1; phi2];
d0 = [d01; d02];
rpar = [rpar1; rpar2];
rperp = [rperp1; rperp2];
w = [w1; w2];

dtod = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; w'];
dtod = [n; dtod(:)];
