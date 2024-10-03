clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/Dropbox';

%ExpNam = 'DiffVACSY'; expno = 128;
%Dguess = [1e-11 2e-9 .8e-10 .8e-10 0.5];
% ExpNam = 'qVAS_Starch'; expno = 65;
% Dguess = [3e-9 3e-9 1e-11 1e-9 0.95];
ExpNam = 'AOTisooctane'; expno = 32; expno = 205;
Dguess = [3e-10 3e-10 2.8633e-11   3.7858e-12   8.9778e-01];
% ExpNam = 'AOTdecanol'; expno = 19;
% Dguess = [3e-10 3e-10 2.8633e-11   3.7858e-12   8.9778e-01];

eval(['load ' DataDir '/' ExpNam '/' num2str(expno) '/NMRacqus'])
eval(['load ' DataDir '/' ExpNam '/' num2str(expno) '/FitDat'])
%figure(1), clf, semilogy(FitDat.Xin{2},FitDat.Yin{2},'o'), return

npeak_water = length(FitDat.Y0);
I = 1*FitDat.Yin{1}/FitDat.Y0(1) + 1*FitDat.Yin{npeak_water}/FitDat.Y0(npeak_water);
I = FitDat.Yin{1}/FitDat.Y0(1);
Iout = 1*FitDat.Yout{1}/FitDat.Y0(1) + 1*FitDat.Yout{npeak_water}/FitDat.Y0(npeak_water);
b = FitDat.Xin{1};
zeta = FitDat.Xin2{1};

Xin = b; Xin2 = zeta;
Yin = I; 
figure(1), clf, semilogy(Xin,Yin,'o'), return

Pin = [Yin(1)*1.05 Dguess]; Funam = 'fDiffVACSYbiexp';

Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = Pin; 
Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,zeros(size(Pin)),[],[],Xin2,Pnorm,Xnorm,Ynorm);
Yout = feval(Funam,Pout,Xin,Xin2); error = Yin - Yout;
%figure(1), clf, semilogy(Xin,Yin,'o',Xin,Yout,'-'), return
%%
FitDat.Pout.Dpar_f = Pout(2);
FitDat.Pout.Dperp_f = Pout(3);
FitDat.Pout.Dpar_s = Pout(4);
FitDat.Pout.Dperp_s = Pout(5);
FitDat.Pout.f = Pout(6);

figure(1), clf
semilogy(Xin,Yin,'o',Xin,Yout,'-')
axis tight

Ib.b = b;
Ib.Deltab = (3*cos(zeta).^2-1)/2;
Ib.I = I;
Ib.Icalc = Yout;

FitDat.Ib = Ib;

eval(['save ' DataDir '/' ExpNam '/' num2str(expno) '/BiexpFitDat FitDat'])

cd(wd)