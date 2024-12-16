function dtor1r2d_nodes_out = dtor1r2d_nodes_select(dtor1r2d_nodes,ind)

[n,par,perp,theta,phi,d0,rpar,rperp,r1,r2] = dtor1r2d_nodes2par(dtor1r2d_nodes);

n = numel(ind);
par = par(ind);
perp = perp(ind);
theta = theta(ind);
phi = phi(ind);
d0 = d0(ind);
rpar = rpar(ind);
rperp = rperp(ind);
r1 = r1(ind);
r2 = r2(ind);

dtor1r2d_nodes_out = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; r1'; r2'];
dtor1r2d_nodes_out = [n; dtor1r2d_nodes_out(:)];
