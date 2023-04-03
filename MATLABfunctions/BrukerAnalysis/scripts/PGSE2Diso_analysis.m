clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
cd(DataDir)

ExpNam = 'isoaniso_test'; expno = 107;

eval(['load ' DataDir '/' ExpNam '/' num2str(expno) '/Dat1D'])
eval(['load ' DataDir '/' ExpNam '/' num2str(expno) '/Dat2D'])

figure(1), clf
semilogy(Dat1D.biso,Dat1D.Iiso,'o',Dat1D.b,Dat1D.I,'o')


b = Dat1D.b; I = Dat1D.I;
%b = Dat1D.biso; I = Dat1D.Iiso;


Xin2 = Dat2D.b;
Xin1 = Dat2D.biso;
Yin = Dat2D.I;
Pin = [1 1e-11 1e-9];
Ycalc = fisoaniso3D(Pin,Xin1,Xin2);

Pin = [1.3 .1e-9 2e-9 1e-10 1e-10 .5];
Funam = 'fisoaniso2comp';

Pout = Pin; Ynorm = mean(reshape(Yin,numel(Yin),1)); Xnorm = mean(reshape(Xin1,numel(Xin1),1)); Pnorm = Pin; 
Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin1/Xnorm,Yin/Ynorm,zeros(size(Pin)),[],[],Xin2/Xnorm,Pnorm,Xnorm,Ynorm);
Yout = feval(Funam,Pout,Xin1,Xin2,ones(size(Pin)),1,1); error = Yin - Yout;
%Ycalc = fisoaniso2comp(Pin,Xin1,Xin2)
chisq = sum(reshape(error.^2,numel(error),1))/numel(error)

figure(2), clf
subplot(2,2,1)
surf(log10(Xin1),log10(Xin2),Yin)
%surf(log10(Xin1),log10(Xin2),real(log10(Ycalc)))
view(0,90), shading('flat'), axis('tight'), axis('square')
subplot(2,2,2)
surf(log10(Xin1),log10(Xin2),Yout)
%surf(log10(Xin1),log10(Xin2),real(log10(Ycalc)))
view(0,90), shading('flat'), axis('tight'), axis('square')
subplot(2,2,3)
surf(log10(Xin1),log10(Xin2),error)
%surf(log10(Xin1),log10(Xin2),real(log10(Ycalc)))
view(0,90), shading('flat'), axis('tight'), axis('square')

figure(1)
hold on
semilogy(Xin1(:,1),Yout(:,1),'-',Xin2(1,:),Yout(1,:),'-')
set(gca,'YLim',[.1 1.2])

NMC = 1000;
if NMC > 0
    noiselevel = std(reshape(error,numel(error),1));
    Pin = Pout; Yin1 = Yin;
    Poutarray = zeros(NMC,length(Pin));
    for nMC = 1:NMC
        nMC
        Yin = Yin1 + noiselevel*randn(size(Yin1));
        Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin1/Xnorm,Yin/Ynorm,zeros(size(Pin)),[],[],Xin2/Xnorm,Pnorm,Xnorm,Ynorm);
        Poutarray(nMC,:) = Pout;
    end
end

%%
figure(3), clf
for nplot = 1:6
    subplot(2,3,nplot)
    hist(Poutarray(:,nplot))
    
    title([num2str(mean(Poutarray(:,nplot)),3) ' ' num2str(std(Poutarray(:,nplot)),1)])
end

%%
MCdat.Poutarray = Poutarray;
MCdat.NMC = NMC;
MCdat.noiselevel = noiselevel;
MCdat.Funam = Funam;

eval(['save ' DataDir '/' ExpNam '/' num2str(expno) '/MCdat MCdat'])


return


mu = .1; Inorm = 1;
lambda = 50;
ND = 50;

%figure(1), clf, semilogy(b,I,'o'), return

D = logspace(-10.5,-8.5,ND)';
bcalc = logspace(7,12,50)';

[barray,Darray] = ndgrid(b,D);
A = exp(-barray.*Darray);
[bcalcarray,Dcalcarray] = ndgrid(bcalc,D);
Acalc = exp(-bcalcarray./Dcalcarray);

L = [1 zeros(1,(ND-1)); zeros(1,(ND-1)) 1];
l = zeros(2,1);

H = [eye((ND-2),ND) zeros((ND-2),2)] + [zeros((ND-2),1) -2*eye((ND-2),ND) zeros((ND-2),1)] + [zeros((ND-2),2) eye((ND-2),ND)];
H(:,(ND+1):(ND+2))=[];
f = zeros((ND-2),1);

PD = Inorm*lsqnonneg([A; mu*H; lambda*L],[I/Inorm; f; l]);
PD(find(isnan(PD)==1)) = eps;
Icalc = A*PD;
error = I - Icalc;

I0 = sum(PD);
D1 = sum(D.*PD)/I0
D2 = sum((D-D1).^2.*PD)/I0

Iinit = I0*exp(-b*D1);

chisq = sum((Icalc-I).^2);
chisqsmooth = sum((f - mu*H*PD).^2);
chisqedge = sum((l - lambda*L*PD).^2);


figure(1), clf
axes('position',[0.1 0.27 0.8 0.65])
semilogy(b,I/I0,'o',b,Icalc/I0,b,Iinit/I0), grid
title(['D1 = ' num2str(D1) ' D2 = ' num2str(D2)])
xlabel('b'), ylabel('I/I_0')
axis([min(b) max(b) min(Icalc/I0) max(I/I0)])
axes('position',[0.58 0.58 0.3 0.3])
semilogx(D,PD/I0,'-ob',[D1 D1],max(PD/I0)*[-.2 1.2],'r')
axis([min(D) max(D) max(PD+eps)/I0*[-.2 1.2]]), grid
xlabel('D'), ylabel('P')
axes('position',[0.1 0.05 0.8 0.1])
plot(b,error,'o'), grid
ylabel('residual')
title(['mean error = ' num2str(mean(error)) '   std error = ' num2str(std(error))])
pause(.1)
