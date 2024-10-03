function Y = fFEXSYfit(Pin,Xin,tmix,Pnorm,Xnorm,Ynorm)

Pin = Pin.*Pnorm;
b = Xin(:,1)'*Xnorm;

Dfast = Pin(1);
Dslow = Pin(2);
Pfast0 = Pin(3);
Pfasteq = Pin(4);
R = Pin(5);
I0vector = Pin(6:length(Pin));

Pfast = [Pfasteq-(Pfasteq - Pfast0)*exp(-R*tmix)]; Pfast(1) = Pfasteq;
%figure(1), clf, semilogx(tmix,Pfast,'-o',tmix,1-Pfast,'-o'), return

Efast = exp(-b*Dfast);
Eslow = exp(-b*Dslow);
%figure(1), clf, semilogx(b,Efast,b,Eslow), return

[Efastarray,Pfastarray] = ndgrid(Efast,Pfast);
[Eslowarray,I0array] = ndgrid(Eslow,I0vector);

Iarray = I0array.*(Pfastarray.*Efastarray + (1-Pfastarray).*Eslowarray);
%figure(1), clf, semilogy(b,Iarray,'-o'), return
Y = Iarray/Ynorm;
