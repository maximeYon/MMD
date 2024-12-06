function Y = fFEXSY4site(Pin,Xin,tmix,bfilt,Pnorm,Xnorm,Ynorm)

Pin = Pin.*Pnorm;
b = Xin(:,1)*Xnorm;

D = Pin(1:4)';
Xeq = Pin(5:7)';
Rvect = Pin(8:13);
I0vector = Pin(14:length(Pin));

Xeq(4,1) = 1 - sum(Xeq);
RMat = [
    0 Rvect(1) Rvect(2) Rvect(3)
    0 0 Rvect(4) Rvect(5)
    0 0 0 Rvect(6)
    0 0 0 0];

t = tmix(2:length(tmix));
Nsites = length(Xeq);
Nbs = length(b);
Ntimes = length(t);

%calculations

X0 = Xeq.*exp(-bfilt*D);
X0 = X0/sum(X0); %X after filter

[bMat,DMat] = ndgrid(b,D);
DiffKernel = exp(-bMat.*DMat);

Ibeq = DiffKernel*Xeq; %equilibrium echo intensity I
Ib0 = DiffKernel*X0; %I after filter, no exchange

[XeqMat1,XeqMat2] = ndgrid(Xeq,Xeq);
KeqMat = XeqMat2./XeqMat1; %matrix of equilibrium constants
%figure(1), clf, surf(XeqMat2.*XeqMat1), view(0,90), shading('flat'), axis('square'), return

[DMat1,DMat2] = ndgrid(D,D);

KexMat = zeros(Nsites,Nsites); %matrix of exchange rate constants
for nsite = 1:Nsites-1
    KexMat(nsite,(nsite+1):Nsites) = RMat(nsite,(nsite+1):Nsites)./(1+KeqMat(nsite,(nsite+1):Nsites));
end
%figure(1), clf, surf(KexMat'), view(0,90), shading('flat'), axis('square'), set(gca,'YDir','normal'), return
KexMat = KexMat + (KexMat.*KeqMat)';
%figure(1), clf, surf(KexMat + KexMat'), view(0,90), shading('flat'), axis('square'), return
KexMat = KexMat - diag(sum(KexMat,1),0);

[EigenVect,EigenValMat] = eig(KexMat); %diagonalize

XtArray = zeros(Nsites,Ntimes);
IbArray = zeros(Nbs,Ntimes);
for ntime = 1:Ntimes
    Xt = (EigenVect*diag(diag(exp(EigenValMat*t(ntime))))*inv(EigenVect))*X0;
    XtArray(:,ntime) = Xt;
    Ib = DiffKernel*Xt; %I for each mixing time
    IbArray(:,ntime) = I0vector(ntime+1)*Ib;
end

%figure(1), clf, semilogy(b,IbArray,'-o'), return

Y = [Ibeq*I0vector(1) IbArray]/Ynorm;
