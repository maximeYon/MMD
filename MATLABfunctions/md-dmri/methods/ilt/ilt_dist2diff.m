function [dtd_nx6,w] = ilt_dist2diff(dtd)

[n,par,perp,theta,phi,w] = ilt_dist2par(dtd);

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
