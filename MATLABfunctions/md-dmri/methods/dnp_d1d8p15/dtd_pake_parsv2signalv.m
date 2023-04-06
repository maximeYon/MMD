function s = dtd_pake_parsv2signalv(dtds, b, bdelta)
% function s = dtd_pake_parsv2signalv(dtds, b, bdelta)
%

a = 3*b.*dtds.diso.*bdelta.*dtds.ddelta + eps;

s = exp(-b.*dtds.diso).*exp(a/3).*...
    sqrt(pi)/2.*real(gammainc(a,1/2)./sqrt(a));

indx = abs(a) < 100*eps;
s(indx) = exp(-b.*dtds.diso(indx));

