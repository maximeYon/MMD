function dtor1r2d_nodes = dtor1r2d_dist2nodes(dtor1r2d)

[n,par,perp,theta,phi,d0,rpar,rperp,r1,r2,w] = dtor1r2d_dist2par(dtor1r2d);

dtor1r2d_nodes = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; r1'; r2'];
dtor1r2d_nodes = [n; dtor1r2d_nodes(:)];
