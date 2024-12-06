function Y = fFEXSYfit(Pin,Xin,tmix,Pnorm,Xnorm,Ynorm)

Pin = Pin.*Pnorm;
b = Xin(:,1)'*Xnorm;

ADC0 = Pin(1);
sigma = Pin(2);
AXR = Pin(3);
I0vector = Pin(4:length(Pin));

[dummy,tmix] = ndgrid(b,tmix);
[b,I0array] = ndgrid(b,I0vector);
ADC = ADC0*(1 - sigma.*exp(-AXR*tmix));
Iarray = exp(-ADC.*b);
Iarray(:,1) = exp(-ADC0.*b(:,1));

Iarray = I0array.*Iarray;

Y = Iarray/Ynorm;
