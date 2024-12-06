function chisq = fchisqFEXSY(Pin)
global Xin Yin weight Pnorm

Pin = Pin.*Pnorm;

Pfasteq = Pin(1);
Pfast0 = Pin(2);
R = Pin(3);

Pfastcalc = Pfasteq-(Pfasteq - Pfast0)*exp(-R*Xin);
Pfastcalc(1) = Pfasteq;

Ycalc = Pfastcalc;
        
chisq = sum((Yin - Ycalc).^2./weight);