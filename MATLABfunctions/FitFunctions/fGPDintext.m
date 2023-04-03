function Y = fGPDKargerfit2(Pin,Xin,deltaarray,tdiffarray,Pnorm,Xnorm,Ynorm)

Pin = Pin.*Pnorm;
qvector = Xin*Xnorm;

D10 = Pin(1);
D2 = Pin(2);
R = Pin(3);
P1 = Pin(4);
%I0vector = Pin(6:length(Pin));

P2 = 1-P1;

gamma = 26.75e7;

alphaR = linspace(0,10,1000);
j = 1./alphaR.*bessel(3/2,alphaR) - bessel(5/2,alphaR);
jj = j(1:length(j)-1).*j(2:length(j));
pos = find(jj<0);
%figure(1), clf, plot(alphaR,j,alphaR(pos),zeros(size(pos)),'o'), grid, return
alphaR = alphaR(pos);

Nqs = length(qvector);
NDds = length(deltaarray);

for nDd = 1:NDds
            delta = deltaarray(nDd);
            tdiff = tdiffarray(nDd);
            Delta = tdiff + delta/3;

            qmax = qvector;
            q2 = 4*pi^2*qmax.^2;

            D1 = fDappSphGPD(delta,Delta,D10,R,alphaR);
 
            E1 = P1.*exp(-q2.*D1.*tdiff);
            E2 = P2.*exp(-q2.*D2.*tdiff);

            E = E1 + E2;

            E1array(:,nDd) = E1;
            E2array(:,nDd) = E2;
end
    
Earray = E1array + E2array;

Y = Earray/Ynorm;

%loglog(qvector,Y,'k-'), pause(.1)
