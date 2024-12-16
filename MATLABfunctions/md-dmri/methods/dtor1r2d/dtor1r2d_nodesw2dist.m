function dtor1r2d = dtor1r2d_nodesw2dist(dtor1r2d_nodes,w)

[n,par,perp,theta,phi,d0,rpar,rperp,r1,r2] = dtor1r2d_nodes2par(dtor1r2d_nodes);

dtor1r2d = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; r1'; r2'; w'];
dtor1r2d = [n; dtor1r2d(:)];
