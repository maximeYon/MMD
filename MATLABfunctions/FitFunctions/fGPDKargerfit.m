function Y = fGPDKargerfit(Pin,Xin,deltaarray,tdiffarray,Pnorm,Xnorm,Ynorm)

Pin = Pin.*Pnorm;
qvector = Xin*Xnorm;

D1 = Pin(1);
D2 = Pin(2);
R = Pin(3);
P1 = Pin(4);
k = Pin(5);
%I0vector = Pin(5:length(Pin));

P2 = 1-P1;
k1 = k/P1;
k2 = k/P2;

dt = 100e-6;
gamma = 26.75e7;

alphaR = linspace(0,10,1000);
j = 1./alphaR.*bessel(3/2,alphaR) - bessel(5/2,alphaR);
jj = j(1:length(j)-1).*j(2:length(j));
pos = find(jj<0);
%figure(1), clf, plot(alphaR,j,alphaR(pos),zeros(size(pos)),'o'), grid, return
alphaR = alphaR(pos);

Nqs = length(qvector);
NDds = length(deltaarray);
Nrep = 1;

for nDd = 1:NDds
            tdiff = tdiffarray(nDd);
            delta = deltaarray(nDd);
            Ndelta = ceil(delta/dt);
            delta = Ndelta*dt;

            Deltaprim = tdiff + delta/3 - delta;
            NDeltaprim = ceil(Deltaprim/dt);
            Deltaprim = NDeltaprim*dt;

            Delta = delta + Deltaprim;
            tdiff = Delta-delta/3;

            D1app = fDappSphGPD(delta,Delta,D1,R,alphaR);
 
            Nt = 2*Ndelta + NDeltaprim + 1;
            t = dt*(0:Nt-1)';

            qmax = qvector;

            gammaGmax = 2*pi*qmax/delta;

            [t, gammaG] = ndgrid(t, gammaGmax);
            
            gammaG(Ndelta+1:Nt-1-Ndelta-1,:) = 0;
            gammaG(Nt-1-Ndelta:Nt-1,:) = -gammaG(Nt-1-Ndelta:Nt-1,:);
            gammaG(Nt,:) = 0;
            q = cumsum(gammaG*dt,1);

            qunit = q;
            tunit = t;

            for n = 1:Nrep-1
                t = [t; tunit + dt + max(t)];
                q = [q; qunit];
            end

            E1 = P1*ones(Nt*Nrep,Nqs);
            E2 = P2*ones(Nt*Nrep,Nqs);

            E1rateDiff = zeros(Nt*Nrep,Nqs);
            E2rateDiff = zeros(Nt*Nrep,Nqs);
            E1rateExch = zeros(Nt*Nrep,Nqs);
            E2rateExch = zeros(Nt*Nrep,Nqs);

            for nt = 2:(Nt*Nrep)
                E1rateDiff(nt,:) = q(nt-1,:).^2*D1app.*E1(nt-1,:);
                E2rateDiff(nt,:) = q(nt-1,:).^2*D2.*E2(nt-1,:);

                E1(nt,:) = E1(nt-1,:) - E1rateDiff(nt,:)*dt;
                E2(nt,:) = E2(nt-1,:) - E2rateDiff(nt,:)*dt;

                E1rateExch(nt,:) = k1*E1(nt,:);
                E2rateExch(nt,:) = k2*E2(nt,:);

                E1(nt,:) = E1(nt,:) - E1rateExch(nt,:)*dt + E2rateExch(nt,:)*dt;
                E2(nt,:) = E2(nt,:) - E2rateExch(nt,:)*dt + E1rateExch(nt,:)*dt;
            end

            E = E1 + E2;

    %         figure(1), clf
    %         subplot(2,1,1)
    %         plot(t,E1rateDiff,t,E2rateDiff)
    %         hold on
    %         plot(t,E1rateExch,'--',t,E2rateExch,'--')
    %         hold off
    %         plot(t(1:Nt*Nrep-1),k)
    %         subplot(2,1,2)
    %         semilogy(t,E1,t,E2,t,E)
    %         title(['E2/E=' num2str(E2(length(E))/E(length(E)),3) ' E=' num2str(E(length(E)),3)])
    %         pause(.1)

            E1array(:,nDd) = E1(Nt,:);
            E2array(:,nDd) = E2(Nt,:);
%             deltavector(nDd) = delta;
%             Deltavector(nDd) = Delta;
%             tdiffvector(nDd) = tdiff;
%             barray(:,nDd) = (gammaGmax*delta).^2*tdiff;
%             qarray(:,nDd) = qmax;
%             Garray(:,nDd) = gammaGmax/gamma;
%             Dapparray(nDd) = D1;
end
    
Earray = E1array + E2array;

Y = Earray/Ynorm;

%loglog(qvector,Y,'k-'), pause(.1)
