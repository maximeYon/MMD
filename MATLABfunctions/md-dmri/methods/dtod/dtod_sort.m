function dtod_out = dtod_sort(dtod_in)

[n,par,perp,theta,phi,d0,rpar,rperp,w] = dtod_dist2par(dtod_in);

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
w = w(ind);

dtod_out = [par'; perp'; theta'; phi'; d0'; rpar'; rperp'; w'];
dtod_out = [n; dtod_out(:)];
