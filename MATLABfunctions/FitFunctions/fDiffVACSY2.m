function Y = fDiffVACSY2(Pin,Xin,Xin2,Pnorm,Xnorm,Ynorm)

if nargin == 3
    Pnorm = ones(size(Pin));
    Xnorm = 1;
    Ynorm = 1;
end

Pin = Pin.*Pnorm;
Xin = Xin*Xnorm;

I0 = Pin(1);
Diso = Pin(2);
DeltaD = Pin(3);

b = Xin;
Deltab = Xin2;

bDDelDel = b.*Diso.*Deltab.*DeltaD;
%max(max(bDDelDel))
%min(min(bDDelDel))

E = exp(-b.*Diso).*exp(bDDelDel).*...
    sqrt(pi)/2.*real(gammainc(3*bDDelDel,1/2)./sqrt(3*bDDelDel));

% E = exp(-b.*Diso).*exp(-2*bDDelDel).*...
%     gammainc(3*bDDelDel,1/2,'scaledlower');
% 
% E = real(E);


indx = bDDelDel == 0;
E(indx) = exp(-b(indx).*Diso);
E(b == 0) = 1;

I = I0*E;

Y = I/Ynorm;

