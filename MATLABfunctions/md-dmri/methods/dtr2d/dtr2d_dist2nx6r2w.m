function [dtd_nx6,r2,w] = dtr2d_dist2nx6r2w(dtr2d)

[n,par,perp,theta,phi,r2,w] = dtr2d_dist2par(dtr2d);

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
