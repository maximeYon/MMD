clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/opt/topspin2/data/DT/nmr';

ExpNam = 'DiffVACSY'; expno = 128;

eval(['load ' DataDir '/' ExpNam '/' num2str(expno) '/NMRacqus'])
eval(['load ' DataDir '/' ExpNam '/' num2str(expno) '/FitDat'])

npeak_water = length(FitDat.Y0);
I = 1*FitDat.Yin{1}/FitDat.Y0(1) + 1*FitDat.Yin{npeak_water}/FitDat.Y0(npeak_water);
Iout = 1*FitDat.Yout{1}/FitDat.Y0(1) + 1*FitDat.Yout{npeak_water}/FitDat.Y0(npeak_water);
b = FitDat.Xin{1};
zeta = FitDat.Xin2{1};

figure(1), clf
semilogy(b,I,'o')
axis tight

PD.minD = 3e-11;
PD.maxD = 5e-9;
PD.NMD = 40;
PD.NDD = 22;

PD.MD = logspace(log10(PD.minD),log10(PD.maxD),PD.NMD);
PD.DD = linspace(-.5,1,PD.NDD);

[PD.MDarDD,PD.DDarMD] = ndgrid(PD.MD,PD.DD);
PD.MDarDDve = reshape(PD.MDarDD,1,PD.NMD*PD.NDD);
PD.DDarMDve = reshape(PD.DDarMD,1,PD.NMD*PD.NDD);

[K.b,K.MD] = ndgrid(b,PD.MDarDDve);
[K.zeta,K.DD] = ndgrid(zeta,PD.DDarMDve);

I0 = 1;

Dpar_zeta = K.MD.*(1 + 2*K.DD.*(3*cos(K.zeta).^2-1)/2);
Dperp_zeta = K.MD.*(1 - K.DD.*(3*cos(K.zeta).^2-1)/2);

E_prolate = zeros(size(K.b));
E_oblate = E_prolate;
E_iso = E_prolate;

E_temp = sqrt(pi)/2*exp(-K.b.*Dperp_zeta)./sqrt(K.b.*(Dpar_zeta-Dperp_zeta)).*erf(sqrt(K.b.*abs(Dpar_zeta-Dperp_zeta)));
indx = Dpar_zeta>Dperp_zeta;
E_prolate(indx) = E_temp(indx);

E_temp = sqrt(pi)/2*exp(-K.b.*Dperp_zeta)./sqrt(K.b.*(Dpar_zeta-Dperp_zeta))*i.*erfi(imag(sqrt(K.b.*(Dpar_zeta-Dperp_zeta))));
indx = Dpar_zeta<Dperp_zeta;
E_oblate(indx) = E_temp(indx);

E_temp = exp(-K.b.*Dperp_zeta);
indx = Dpar_zeta==Dperp_zeta;
E_iso(indx) = E_temp(indx);

E = E_prolate + E_oblate + E_iso;

A = I0*E;

figure(2), clf
semilogy(K.b,A,'-')
set(gca,'YLim',[1e-2 1.1])


lambda = 1e-3;
lambda = 1e-1;
H = 2*A.'*A;
f = lambda*ones(PD.NMD*PD.NDD,1)-2*(I.'*A).';

lb = zeros(PD.NMD*PD.NDD,1);
%lb = [];

PQP = quadprog(H,f,[],[],[],[],lb,[]);
%PQP = quadprog(H,f,-eye(PD.NMD*PD.NDD),eps*ones(PD.NMD*PD.NDD,1),[],[],[],[]);
%options = optimset('MaxIter',1e4,'LargeScale','off');
%PQP = lsqlin(A,I,-eye(PD.NMD*PD.NDD),eps*ones(PD.NMD*PD.NDD,1),[],[],[],[],[],options);
%PQP = lsqlin(A,I,[],[],[],[],lb,[],[],options);


Ifit = A*PQP;
epsilon = (Ifit - I).^2;
chisq = sum(epsilon)
l1 = lambda*sum(abs(PQP))
TV = sum(abs(gradient(PQP)))

figure(3), clf
semilogy(b,I,'o',b,Ifit,'-',b,Iout,'-')
axis tight

epsilon2 = (Iout - I).^2;
chisq2 = sum(epsilon2)

PD.I = reshape(PQP,PD.NMD,PD.NDD);
%%
figure(4), clf
imagesc(log10(PD.MD),PD.DD,PD.I')
set(gca,'YDir','normal')
% surf(log10(PD.MDarDD),PD.DDarMD,PD.I)
% view(30,90)
shading interp
axis tight
colormap('jet')

Ib.b = b;
Ib.Deltab = (3*cos(zeta).^2-1)/2;
Ib.I = I;

eval(['save ' DataDir '/' ExpNam '/' num2str(expno) '/ILTdat PD Ib'])

cd(wd)