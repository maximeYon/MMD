function Y = fGPDKargerfit2(Pin,Xin,deltaarray,tdiffarray,Pnorm,Xnorm,Ynorm)

Pin = Pin.*Pnorm;
qvector = Xin*Xnorm;

D10 = Pin(1);
D2 = Pin(2);
R = Pin(3);
P1 = Pin(4);
H = Pin(5);
%I0vector = Pin(6:length(Pin));

P2 = 1-P1;
k1 = 3*H/R;
k2 = k1*P1/P2;

gamma = 26.75e7;

alphaR = linspace(0,10,1000);
%j = 1./alphaR.*bessel(3/2,alphaR) - bessel(5/2,alphaR);
j = 1./alphaR.*besselj(3/2,alphaR) - besselj(5/2,alphaR);
jj = j(1:length(j)-1).*j(2:length(j));
pos = find(jj<0);
%figure(1), clf, plot(alphaR,j,alphaR(pos),zeros(size(pos)),'o'), grid, return
alphaR = alphaR(pos);

Nqs = length(qvector);
NDds = length(deltaarray);
Nrep = 1;

for nDd = 1:NDds
            delta = deltaarray(nDd);
            tdiff = tdiffarray(nDd);
            Delta = tdiff + delta/3;

            D1 = fDappSphGPD(delta,Delta,D10,R,alphaR);
 
            qmax = qvector;

            gammaGmax = 2*pi*qmax/delta;

            qmax = qvector;
            q2 = 4*pi^2*qmax.^2;

            gammaGmax = 2*pi*qmax/delta;

            DA = 0.5*(D1+D2+1./q2*(k1+k2) - sqrt((D1-D2+1./q2*(k1-k2)).^2 + 4./q2.^2*k1*k2));
            DB = 0.5*(D1+D2+1./q2*(k1+k2) + sqrt((D1-D2+1./q2*(k1-k2)).^2 + 4./q2.^2*k1*k2));

            PB = (P1*D1+P2*D2-DA)./(DB-DA);
            PA = 1-PB;

            E1 = PA.*exp(-q2.*DA.*tdiff);
            E2 = PB.*exp(-q2.*DB.*tdiff);

            E = E1 + E2;

            E1array(:,nDd) = E1;
            E2array(:,nDd) = E2;
end
    
Earray = E1array + E2array;

Y = Earray/Ynorm;

%loglog(qvector,Y,'k-'), pause(.1)
