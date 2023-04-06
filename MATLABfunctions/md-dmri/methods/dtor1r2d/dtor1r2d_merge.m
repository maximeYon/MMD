function dtor1r2d = dtor1r2d_merge(dtor1r2d1,dtor1r2d2)

[n1,par1,perp1,theta1,phi1,d01,rpar1,rperp1,r11,r21,w1] = dtor1r2d_dist2par(dtor1r2d1);
[n2,par2,perp2,theta2,phi2,d02,rpar2,rperp2,r12,r22,w2] = dtor1r2d_dist2par(dtor1r2d2);

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
w = [w1; w2];

dtor1r2d = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; r1'; r2'; w'];
dtor1r2d = [n; dtor1r2d(:)];
