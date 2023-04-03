function Y = fGPDintext2(Pin,Xin,deltaarray,tdiffarray,Pnorm,Xnorm,Ynorm)

Pin = Pin.*Pnorm;
qvector = Xin*Xnorm;

D2 = Pin(1);
P1 = Pin(2);
D1vector = Pin(3:length(Pin));
%I0vector = Pin(6:length(Pin));

P2 = 1-P1;

NDds = length(deltaarray);

for nDd = 1:NDds
            delta = deltaarray(nDd);
            tdiff = tdiffarray(nDd);
            Delta = tdiff + delta/3;

            qmax = qvector;
            q2 = 4*pi^2*qmax.^2;

            D1 = D1vector(nDd);
 
            E1 = P1.*exp(-q2.*D1.*tdiff);
            E2 = P2.*exp(-q2.*D2.*tdiff);

            E = E1 + E2;

            E1array(:,nDd) = E1;
            E2array(:,nDd) = E2;
end
    
Earray = E1array + E2array;

Y = Earray/Ynorm;

%loglog(qvector,Y,'k-'), pause(.1)
