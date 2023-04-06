function dtor1r2d_out = dtor1r2d_sort(dtor1r2d_in)

[n,par,perp,theta,phi,d0,rpar,rperp,r1,r2,w] = dtor1r2d_dist2par(dtor1r2d_in);

[wsort,ind] = sort(w,'descend');
ind = ind(wsort>0);

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
w = w(ind);

dtor1r2d_out = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; r1'; r2'; w'];
dtor1r2d_out = [n; dtor1r2d_out(:)];
