function Y = fFEXSYfitAXR3(Pin,b,Para,Pnorm,Xnorm,Ynorm)

Pin = Pin.*Pnorm;
b = b(:,1)*Xnorm;
tm = Para.tm;
bf = Para.bf;
Nb = length(b);

D1 = Pin(1);
D2 = Pin(2);
f1eq = Pin(3);
AXR = Pin(4);
I0 = Pin(5);
I1tmvector = [1 1 Pin(6:length(Pin))];

[b,tm] = ndgrid(b,tm);

f2eq = 1-f1eq;
ADC = f1eq*D1 + f2eq*D2;
reldiffD = (D1-D2)/ADC;

f10 = f1eq*exp(-bf*D1)/(f1eq*exp(-bf*D1)+f2eq*exp(-bf*D2));
f20 = 1-f10;

sigma = reldiffD*(f1eq-f10);

ADCtm = ADC*(1-sigma*exp(-AXR*tm));

I1 = I0*(f1eq*exp(-bf*D1) + f2eq*exp(-bf*D2));
I1tm = I1*I1tmvector;
[dummy,I1tm] = ndgrid(1:Nb,I1tm);

ADCtm(:,1) = ADC*ones(Nb,1);
I1tm(:,1) = I0*ones(Nb,1);

Iarray = I1tm.*exp(-b.*ADCtm);

Y = Iarray/Ynorm;
