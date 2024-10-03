function Y = fFEXSYfit(Pin,Xin,tmix,Pnorm,Xnorm,Ynorm)

Pin = Pin.*Pnorm;
b = Xin(:,1)'*Xnorm;

Dfast = Pin(1);
Dslow = Pin(2);
I0fast = Pin((1:length(tmix)) + 2);
I0slow = Pin((1:length(tmix)) + 2 + length(tmix));

Efast = exp(-b*Dfast);
Eslow = exp(-b*Dslow);
%figure(1), clf, semilogx(b,Efast,b,Eslow), return

[Efastarray,I0fastarray] = ndgrid(Efast,I0fast);
[Eslowarray,I0slowarray] = ndgrid(Eslow,I0slow);

Iarray = I0fastarray.*Efastarray + I0slowarray.*Eslowarray;
%figure(1), clf, semilogy(b,Iarray,'-o'), return
Y = Iarray/Ynorm;
