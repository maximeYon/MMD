function dtd_out = dtd_pa_sort(dtd_in)

[n,par,perp,w] = dtd_pa_dist2par(dtd_in);

[wsort,ind] = sort(w,'descend');
ind = ind(wsort>0);

n = numel(ind);
par = par(ind);
perp = perp(ind);
w = w(ind);

dtd_out = [par'; perp'; w'];
dtd_out = [n; dtd_out(:)];
