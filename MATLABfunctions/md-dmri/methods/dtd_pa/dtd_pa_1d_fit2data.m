function s = dtd_pa_1d_fit2data(m, xps)
% function s = dtd_pa_1d_fit2data(m, xps)

if (m(1) > 0)
    
    dtd_pa = dtd_pa_m2dtd(m);
    [~,par,perp,w] = dtd_pa_dist2par(dtd_pa);
    diso = (par + 2*perp)/3;
    ddelta = (par - perp)/3./diso;
    
    bd = xps.b*diso';
    bddeldel = (xps.b.*xps.b_delta)*(diso'.*ddelta');
    k = exp(-bd).*exp(bddeldel).*...
        sqrt(pi)/2.*real(gammainc(3*bddeldel,1/2)./sqrt(3*bddeldel));

    % Correct an issue where gammainc(x)/sqrt(x) is undefined for small x
    indx = abs(bddeldel) < 10 * eps;
    k(indx) = exp(-bd(indx));
    k(bd == 0) = 1;
    k(bddeldel < -10) = 0;
    k(isnan(k)) = 0;
    k(isinf(k)) = 0;
    
    s = k*w;
else
    s = zeros(xps.n,1);
end

