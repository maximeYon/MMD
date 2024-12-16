function [dtor1r2d_nx6,r1,r2,w] = dtor1r2d_dist2nx6r1r2w(dtor1r2d,xps)

[n,dpar,dperp,theta,phi,d0,rpar,rperp,r1,r2,w] = dtor1r2d_dist2par(dtor1r2d);
dtor1r2d_nodes = dtor1r2d_dist2nodes(dtor1r2d);
[dtor1r2d_nx6,r1,r2] = dtor1r2d_nodes2nx6r1r2(dtor1r2d_nodes,xps);


return
[n,par,perp,theta,phi,r1,r2,w] = dtor1r2d_dist2par(dtor1r2d);

xcos = cos(phi).*sin(theta);
ycos = sin(phi).*sin(theta);
zcos = cos(theta);

trace = par + 2*perp;
delta = (par - perp)./trace;

xx = trace/3.*(1 + delta.*(3*xcos.*xcos - 1));
xy = trace/3.*(0 + delta.*(3*xcos.*ycos - 0));
xz = trace/3.*(0 + delta.*(3*xcos.*zcos - 0));
yy = trace/3.*(1 + delta.*(3*ycos.*ycos - 1));
yz = trace/3.*(0 + delta.*(3*ycos.*zcos - 0));
zz = trace/3.*(1 + delta.*(3*zcos.*zcos - 1));

dtd_nx6 = [xx yy zz sqrt(2)*[xy xz yz]];
