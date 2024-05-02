function dtd_out = dtd_sort(dtd_in)

[n,par,perp,theta,phi,w] = dtd_dist2par(dtd_in);

[wsort,ind] = sort(w,'descend');
ind = ind(wsort>0);

n = numel(ind);
par = par(ind);
perp = perp(ind);
theta = theta(ind);
phi = phi(ind);
w = w(ind);

dtd_out = [par'; perp'; theta'; phi'; w'];
dtd_out = [n; dtd_out(:)];
