function Y = fFEXSYfitR1(Pin,b,tmix,Pnorm,Xnorm,Ynorm)

Pin = Pin.*Pnorm;
b = b(:,1)*Xnorm;

ADC = Pin(1);
R1 = Pin(2);
A = Pin(3);
I0 = Pin(4);

[b,tmix] = ndgrid(b,tmix);
I0array = I0*(1 - A.*exp(-R1*tmix));
Iarray = I0array.*exp(-ADC.*b);

Y = Iarray/Ynorm;
