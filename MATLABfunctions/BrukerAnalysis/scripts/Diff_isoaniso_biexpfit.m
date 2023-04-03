clear all

wd = cd;

%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/opt/topspin2/data/DT/nmr';
DataDir = '/Users/daniel/Dropbox';

% ExpNam = 'AOTdecanol'; expno = 19;
% Dguess = [1.7983e-11   1.8141e-09   8.2607e-11   8.2810e-11   8.9363e-01];
% %Dguess = [2*1.8066e-09   5.3400e-11   8.3012e-11   8.2962e-11   6.2430e-01];

ExpNam = 'Yeast_tripgse'; expno = 53;
Dguess = [1.8141e-09 1.7983e-11 2e-11   2e-11   .9];
Dguess = [1.8141e-09 .95 2e-11 .1  .9];

eval(['load ' DataDir '/' ExpNam '/' num2str(expno) '/NMRacqus'])
eval(['load ' DataDir '/' ExpNam '/' num2str(expno) '/FitDat'])
%figure(1), clf, semilogy(FitDat.Xin{2},FitDat.Yin{2},'o'), return

npeak = length(FitDat.Y0);
Nb = AcqDat.Nb;

% Xval = FitDat.Xin{npeak};
% X2val = (3*cos(FitDat.Xin2{npeak}).^2 - 1)/2;
% Yinval = (.2*FitDat.Yin{1} + 1*FitDat.Yin{npeak})./(.2*FitDat.Y0(1) + 1*FitDat.Y0(npeak));
% 
% Xval = reshape(Xval,sqrt(Nb),sqrt(Nb));
% X2val = reshape(X2val,sqrt(Nb),sqrt(Nb));
% Yinval = reshape(Yinval,sqrt(Nb),sqrt(Nb));
% 
% bzz = (2*X2val+1)/3.*Xval;
% bxx = (Xval - bzz)/2;
% 
% bmin_array = zeros(sqrt(Nb),sqrt(Nb),2);
% bmin_array(:,:,1) = bzz;
% bmin_array(:,:,2) = bxx;
% 
% biso = 3*min(bmin_array,[],3);
% 
% beA = bzz - biso/3;
% beR = 2*(bxx - biso/3);
% baniso = beA - beR;
% 
% Xval = biso;
% X2val = baniso;
% 
% figure(1), clf
% %subplot(2,1,1)
% %semilogy(Xval,Yinval,'b-')
% semilogy(Xval + X2val,Yinval,'b-')
% %hold on
% %subplot(2,1,2)
% %semilogy((X2val)',Yinval','r-')
% %semilogy((Xval + X2val)',Yinval','r-')
% axis([1.5e10*[-.1 1.1] 1e-1 1.1])

%figure(1), clf, semilogy(FitDat.Xin{2},I,'o'), return

I = FitDat.Yin{npeak}./FitDat.Y0(npeak);
%I = (.2*FitDat.Yin{1} + 1*FitDat.Yin{npeak})./(.2*FitDat.Y0(1) + 1*FitDat.Y0(npeak));
%Iout = (1*FitDat.Yout{1} + 1*FitDat.Yout{npeak})./(1*FitDat.Y0(1) + 1*FitDat.Y0(npeak));
b = FitDat.Xin{1};
zeta = FitDat.Xin2{1};

Xin = b; Xin2 = (3*cos(zeta).^2-1)/2;
Yin = I; 
%figure(1), clf, semilogy(Xin,Yin,'o'), return

Pin = [Yin(1)*1.1 Dguess]; Funam = 'fDiffVACSYbiexp2';
vlb = [Yin(1)*.1 1e-13 -.5 1e-13 -.5 0];
vub = [Yin(1)*1.1 1e-8 1 1e-8 1 1];

Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = abs(Pin); 
Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,vlb./Pnorm,vub./Pnorm,[],Xin2,Pnorm,Xnorm,Ynorm);
Yout = feval(Funam,Pout,Xin,Xin2); error = Yin - Yout;
%figure(1), clf, semilogy(Xin,Yin,'o',Xin,Yout,'+'), return
%%
FitDat.Pout.I0 = Pout(1);
FitDat.Pout.Diso_f = Pout(2);
FitDat.Pout.DeltaD_f = Pout(3);
FitDat.Pout.Diso_s = Pout(4);
FitDat.Pout.DeltaD_s = Pout(5);
FitDat.Pout.Dpar_f = Pout(2)*(1+2*Pout(3));
FitDat.Pout.Dperp_f = Pout(2)*(1-Pout(3));
FitDat.Pout.Dpar_s = Pout(4)*(1+2*Pout(5));
FitDat.Pout.Dperp_s = Pout(4)*(1-Pout(5));
FitDat.Pout.f = Pout(6);

% figure(1), clf
% semilogy(Xin,Yin,'o',Xin,Yout,'+')
% axis tight

Ib.b = b;
Ib.Deltab = (3*cos(zeta).^2-1)/2;
Ib.bzz = (2*Ib.Deltab+1)/3.*Ib.b;
Ib.bxx = (Ib.b - Ib.bzz)/2;

Ib.biso = 3*min([Ib.bzz Ib.bxx],[],2);

Ib.beA = Ib.bzz - Ib.biso/3;
Ib.beR = 2*(Ib.bxx - Ib.biso/3);
Ib.baniso = Ib.beA - Ib.beR;

Ib.I = I;
Ib.Icalc = Yout;

FitDat.Ib = Ib;

bcalc = max(Ib.b)*linspace(0,1,1000)';
Deltabcalc_iso = zeros(size(bcalc));
zetacalc_iso = acos(1/sqrt(3))*ones(size(bcalc));
Deltabcalc_aniso = (bcalc-min(Ib.biso))./bcalc;
zetacalc_aniso = real(acos(sqrt((2*Deltabcalc_aniso+1)/3)));
Deltabcalc_aniso2 = bcalc./(bcalc+min(Ib.biso));
zetacalc_aniso2 = real(acos(sqrt((2*Deltabcalc_aniso2+1)/3)));

Pout = [1 FitDat.Pout.Diso_f FitDat.Pout.DeltaD_f FitDat.Pout.Diso_s FitDat.Pout.DeltaD_s FitDat.Pout.f];
Icalc_iso = fDiffVACSYbiexp2(Pout,bcalc,Deltabcalc_iso);
Icalc_aniso = fDiffVACSYbiexp2(Pout,bcalc,Deltabcalc_aniso);
Icalc_aniso2 = fDiffVACSYbiexp2(Pout,bcalc,Deltabcalc_aniso2);

Xval = 1e-9*reshape(Ib.biso,sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
Yval = 1e-9*reshape(Ib.baniso,sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
Zval = reshape(Ib.I/FitDat.Pout.I0,sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
Zval2 = reshape(Ib.Icalc/FitDat.Pout.I0,sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));

Imin = .008; Imin = 3e-2;
bmax = 25e9; bmax = max(max(Ib.biso));

Zval(Zval<Imin) = NaN;
Zval2(Zval2<Imin) = NaN;
Icalc_aniso2(Icalc_aniso2<Imin) = NaN;
Icalc_iso(Icalc_iso<Imin) = NaN;

lw = 3;
fs = 20;
bottom = .15;
height = .84;


figure(2), clf
axh1 = axes('position',[0.12 bottom .32 height]);
lh12 = semilogy(bcalc*1e-9,Icalc_aniso,'g-');
hold on
lh11 = semilogy(bcalc*1e-9,Icalc_iso,'b-');
lh14 = semilogy((Xval(1,:) + Yval(1,:))',Zval(1,:)','go');
lh13 = semilogy(Xval(:,1),Zval(:,1),'bo');
%lh14 = semilogy((Xval(1,:) + Yval(1,:))'',Zval2(1,:)','k-');
set(axh1,'XLim',bmax*1e-9*[-.1 1.1],'YLim',[Imin 1.1],'TickLength',.02*[1 1])

axh2 = axes('position',[0.58 bottom .4 height]);
lh21 = plot3(Xval,Yval,log10(Zval2),'k-');
hold on
lh22 = plot3(Xval',Yval',log10(Zval2)','k-');
lh23 = plot3(Xval,Yval,log10(Zval),'ko');
lh24 = plot3(bcalc*1e-9,zeros(size(bcalc)),log10(Icalc_iso),'b-');
lh25 = plot3(1e-9*min(Ib.biso)*ones(size(bcalc)),bcalc*1e-9,log10(Icalc_aniso2),'g-');
lh26 = plot3(Xval(:,1),zeros(size(Xval(:,1))),log10(Zval(:,1)),'bo');
lh27 = plot3(Xval(1,:)',Yval(1,:)',log10(Zval(1,:))','go');

view(160,10)
axis([bmax*1e-9*[-.1 1.1 -.1 1.1] log10(Imin) .1])

set([axh1 axh2],'LineWidth',lw,'FontSize',fs,'TickDir','out','Box','off')
Itick = [.01 .1 1]; ItickLabel = {'0.01', '0.1', '1'};
set(axh1,'YTick',Itick,'YTickLabel',ItickLabel)
set(axh2,'ZTick',log10(Itick),'ZTickLabel',ItickLabel)
set([lh11 lh24],'LineWidth',1*lw)
set([lh12 lh25],'LineWidth',1.5*lw,'Color',[0 .7 0])
set([lh13 lh26],'LineWidth',.1*lw,'MarkerSize',3*lw,'MarkerFaceColor',[0 0 1])
set([lh14 lh27],'LineWidth',.1*lw,'MarkerSize',4*lw,'Color',[0 .7 0],'MarkerFaceColor',[0 .7 0])
set([lh21 lh22],'LineWidth',.25*lw)
set([lh23],'LineWidth',.1*lw,'MarkerSize',2*lw,'Color',[0 0 0],'MarkerFaceColor',[0 0 0])

eval(['save ' DataDir '/' ExpNam '/' num2str(expno) '/BiexpFitDat FitDat'])
eval(['print ' DataDir '/' ExpNam '/' num2str(expno) '/BiexpFitPlot -depsc -loose'])

cd(wd)