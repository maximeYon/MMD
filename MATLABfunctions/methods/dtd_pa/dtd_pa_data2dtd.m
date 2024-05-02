function dtd = dtd_pa_data2dtd(stemp,b,b_delta,dtd_nodes)

[n,par,perp] = dtd_pa_nodes2par(dtd_nodes);

diso = (par + 2*perp)/3;
ddelta = (par - perp)/3./diso;

bd = b*diso';
bddeldel = (b.*b_delta)*(diso'.*ddelta');
k = exp(-bd).*exp(bddeldel).*...
    sqrt(pi)/2.*real(gammainc(3*bddeldel,1/2)./sqrt(3*bddeldel));

indx = bddeldel == 0;
k(indx) = exp(-bd(indx));
k(bd == 0) = 1;
k(bddeldel < -10) = 0;
k(isnan(k)) = 0;
k(isinf(k)) = 0;

% dtd_nx6 = dtd_nodes2nx6(dtd_nodes);
% 
% k = exp(-bt_mx6*dtd_nx6');
snorm = max(stemp);
w = snorm*lsqnonneg(k,stemp/snorm);

dtd = dtd_pa_nodesw2dist(dtd_nodes,w);
dtd = dtd_pa_sort(dtd);




