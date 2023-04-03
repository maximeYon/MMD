clear all

fs = 20;

DataDir = '/Users/daniel/Dropbox';
ExpNam = 'AOToct5';

expno = 13;
eval(['load ' DataDir '/' ExpNam '/' num2str(expno) '/Images'])
ImagesR2 = Images;
expno = 11;
eval(['load ' DataDir '/' ExpNam '/' num2str(expno) '/Images'])
ImagesDT = Images;
Nvoxels = numel(ImagesR2.S0);
Nx = length(r.x);
Ny = length(r.y);
Nz = length(r.z);

%%
histdat.R2S0 = reshape(ImagesR2.S0,Nvoxels,1);
histdat.R2 = reshape(ImagesR2.R2,Nvoxels,1);
histdat.DTS0 = reshape(ImagesDT.S0,Nvoxels,1);
histdat.ADC = reshape(ImagesDT.ADC,Nvoxels,1);
histdat.FA = reshape(ImagesDT.FA,Nvoxels,1);
histdat.CL = reshape(ImagesDT.CL,Nvoxels,1);
histdat.CP = reshape(ImagesDT.CP,Nvoxels,1);
histdat.v1.x = reshape(ImagesDT.v1.x,Nvoxels,1);
histdat.v1.y = reshape(ImagesDT.v1.y,Nvoxels,1);
histdat.v1.z = reshape(ImagesDT.v1.z,Nvoxels,1);
histdat.P2 = (3*histdat.v1.z.^2 - 1)/2;
histdat.theta = acos(sqrt((2*histdat.P2 + 1)/3));
histdat.z = reshape(repmat(r.z, [1 Ny Nx]),Nvoxels,1);

mask.S0 = find(histdat.R2S0>0 & histdat.DTS0>0);
mask.S0R2 = find(histdat.R2S0>0 & histdat.DTS0>0 & histdat.R2<300);
mask.S0R2CL = find(histdat.R2S0>0 & histdat.DTS0>0 & histdat.R2<400 & histdat.R2>1 ...
    & histdat.CL>0.7 & histdat.z > -7e-3 & histdat.z < 7e-3 ...
    & histdat.theta > 0/180*pi & histdat.theta < 90/180*pi);

Nsubx = 4;
Nsuby = 2;
Nbins = 500;
figure(1), clf
subplot(Nsuby,Nsubx,1)
hist(histdat.R2(mask.S0R2CL),Nbins)
xlabel('R2')
subplot(Nsuby,Nsubx,2)
hist(histdat.ADC(mask.S0R2CL)*1e9,Nbins)
xlabel('ADC')
subplot(Nsuby,Nsubx,3)
hist(histdat.FA(mask.S0R2CL),Nbins)
xlabel('FA')
subplot(Nsuby,Nsubx,4)
hist(histdat.CL(mask.S0R2CL),Nbins)
xlabel('CL')
subplot(Nsuby,Nsubx,5)
hist(histdat.CP(mask.S0R2CL),Nbins)
xlabel('CP')
subplot(Nsuby,Nsubx,6)
hist(histdat.theta(mask.S0R2CL)/pi*180,Nbins)
xlabel('\theta')

%subplot(Nsuby,2,4)
figure(11), clf
x = histdat.theta(mask.S0R2CL)/pi*180;
[x,indx] = sort(x,1,'ascend');
y = histdat.R2(mask.S0R2CL);
y = y(indx);
p = polyfit(x,y,10);
xfit = linspace(0,90,100)';
yfit = polyval(p,xfit);
yfit = 30 + 150*abs((3*cos(xfit/180*pi).^2-1)/2).^2;
plot(x,y,'.','MarkerSize',5)
hold on, plot(xfit,yfit,'k-','LineWidth',10)
set(gca,'FontSize',50,'LineWidth',5,'Box','off','TickDir','out','TickLength',.03*[1 1])
axis([0 90 0 200])
set(gca,'XTick',[0 30 54.7 90])
%xlabel('\theta'), ylabel('R2')
eval(['print ' DataDir '/' ExpNam '/' num2str(expno) '/CorrR2theta -djpeg -r150 -loose'])

%%

Nslices = 3^2;
slices = round(Nz*linspace(.3,.7,Nslices+2)); slices([1 Nslices+2]) = [];
width = .5/sqrt(Nslices);
height = 11/8.5*.5/sqrt(Nslices);

figure(4), clf           

c.bright = zeros(Nz,Ny,Nx);
c.bright(mask.S0R2CL) = ImagesDT.CL(mask.S0R2CL);
c.r = abs(ImagesDT.v1.x);
c.g = abs(ImagesDT.v1.y);
c.b = abs(ImagesDT.v1.z);

c.bright(isnan(c.bright)) = 0; c.bright(isinf(c.bright)) = 0;
c.r(isnan(c.r)) = 0; c.r(isinf(c.r)) = 0;
c.g(isnan(c.g)) = 0; c.g(isinf(c.g)) = 0;
c.b(isnan(c.b)) = 0; c.b(isinf(c.b)) = 0;
c.bright(c.bright<0) = 0; c.bright(c.bright>1) = 1;
c.r(c.r<0) = 0; c.r(c.r>1) = 1;
c.b(c.b<0) = 0; c.b(c.b>1) = 1;
c.g(c.g<0) = 0; c.g(c.g>1) = 1;

ny = Ny/2 + 1;

axes('position',[.07 .1 .2 .9],'FontSize',fs)
C = zeros(length(r.z),length(r.x),3);
C(:,:,1) = squeeze(c.r(:,ny,:).*c.bright(:,ny,:));
C(:,:,2) = squeeze(c.g(:,ny,:).*c.bright(:,ny,:));
C(:,:,3) = squeeze(c.b(:,ny,:).*c.bright(:,ny,:));
image([min(r.x) max(r.x)]*1e3,[min(r.z) max(r.z)]*1e3,C)
set(gca,'YDir','normal')
axis equal, axis tight
%ylabel('\itz\rm / mm'), xlabel('\itx\rm / mm')

nx = Nx/2 + 1;
axes('position',[.28 .1 .2 .9],'FontSize',fs)
C = zeros(length(r.z),length(r.y),3);
C(:,:,1) = squeeze(c.r(:,:,nx).*c.bright(:,:,nx));
C(:,:,2) = squeeze(c.g(:,:,nx).*c.bright(:,:,nx));
C(:,:,3) = squeeze(c.b(:,:,nx).*c.bright(:,:,nx));
image([min(r.y) max(r.y)]*1e3,[min(r.z) max(r.z)]*1e3,C)
set(gca,'YDir','normal','YTick',[])
axis equal, axis tight
%xlabel('\ity\rm / mm')

nslice = 1;

for county = 1:sqrt(Nslices)
    for countx = 1:sqrt(Nslices)
        axes('position',[1-countx*width .2+(county-1)*height 1.01*width 1.01*height])
        nz = slices(nslice);
        C = zeros(length(r.y),length(r.x),3);
        C(:,:,1) = squeeze(c.r(nz,:,:).*c.bright(nz,:,:));
        C(:,:,2) = squeeze(c.g(nz,:,:).*c.bright(nz,:,:));
        C(:,:,3) = squeeze(c.b(nz,:,:).*c.bright(nz,:,:));
        image([min(r.x) max(r.x)]*1e3,[min(r.y) max(r.y)]*1e3,C)
        set(gca,'YDir','normal','YTick',[])
        axis tight, axis off
        nslice = nslice+1;
    end
end

eval(['print ' DataDir '/' ExpNam '/' num2str(expno) '/CL -depsc -loose'])

figure(5), clf           
%Imax = max(reshape(ImagesR2.S0,Nvoxels,1));
Imax = 6e5;
R2max = 180; R2min = 30;

c.bright = zeros(Nz,Ny,Nx);
c.bright(mask.S0R2CL) = ImagesR2.S0(mask.S0R2CL)/Imax;
C = (log10(ImagesR2.R2) - log10(R2min))/(log10(R2max) - log10(R2min));
%C = (ImagesR2.R2 - R2min)/(R2max - R2min);
C(C<0) = 0; C(C>1) = 1;
c.r = 2*C;
c.g = ones(size(C)) - abs(1-2*C);
c.b = 1-2*C;

c.bright(isnan(c.bright)) = 0; c.bright(isinf(c.bright)) = 0;
c.r(isnan(c.r)) = 0; c.r(isinf(c.r)) = 0;
c.g(isnan(c.g)) = 0; c.g(isinf(c.g)) = 0;
c.b(isnan(c.b)) = 0; c.b(isinf(c.b)) = 0;
c.bright(c.bright<0) = 0; c.bright(c.bright>1) = 1;
c.r(c.r<0) = 0; c.r(c.r>1) = 1;
c.b(c.b<0) = 0; c.b(c.b>1) = 1;
c.g(c.g<0) = 0; c.g(c.g>1) = 1;

ny = Ny/2 + 1;

axes('position',[.07 .1 .2 .9],'FontSize',fs)
C = zeros(length(r.z),length(r.x),3);
C(:,:,1) = squeeze(c.r(:,ny,:).*c.bright(:,ny,:));
C(:,:,2) = squeeze(c.g(:,ny,:).*c.bright(:,ny,:));
C(:,:,3) = squeeze(c.b(:,ny,:).*c.bright(:,ny,:));
image([min(r.x) max(r.x)]*1e3,[min(r.z) max(r.z)]*1e3,C)
set(gca,'YDir','normal')
axis equal, axis tight
%ylabel('\itz\rm / mm'), xlabel('\itx\rm / mm')

nx = Nx/2 + 1;
axes('position',[.28 .1 .2 .9],'FontSize',fs)
C = zeros(length(r.z),length(r.y),3);
C(:,:,1) = squeeze(c.r(:,:,nx).*c.bright(:,:,nx));
C(:,:,2) = squeeze(c.g(:,:,nx).*c.bright(:,:,nx));
C(:,:,3) = squeeze(c.b(:,:,nx).*c.bright(:,:,nx));
image([min(r.y) max(r.y)]*1e3,[min(r.z) max(r.z)]*1e3,C)
set(gca,'YDir','normal','YTick',[])
axis equal, axis tight
%xlabel('\ity\rm / mm')

width = .5/sqrt(Nslices);
height = 11/8.5*.5/sqrt(Nslices);

nslice = 1;

for county = 1:sqrt(Nslices)
    for countx = 1:sqrt(Nslices)
        axes('position',[1-countx*width .2+(county-1)*height 1.01*width 1.01*height])
        nz = slices(nslice);
        C = zeros(length(r.y),length(r.x),3);
        C(:,:,1) = squeeze(c.r(nz,:,:).*c.bright(nz,:,:));
        C(:,:,2) = squeeze(c.g(nz,:,:).*c.bright(nz,:,:));
        C(:,:,3) = squeeze(c.b(nz,:,:).*c.bright(nz,:,:));
        image([min(r.x) max(r.x)]*1e3,[min(r.y) max(r.y)]*1e3,C)
        set(gca,'YDir','normal','YTick',[])
        axis tight, axis off
        nslice = nslice+1;
    end
end

eval(['print ' DataDir '/' ExpNam '/' num2str(expno) '/R2 -depsc -loose'])

%%
%Smooth ODF nodes

load UDSRTriN1000
ODFsmooth = UDSR;

%Discrete ODF
ODFdiscrete.x = histdat.v1.x(mask.S0R2CL);
ODFdiscrete.y = histdat.v1.y(mask.S0R2CL);
ODFdiscrete.z = histdat.v1.z(mask.S0R2CL);
RandSample = ceil(rand(10000,1)*numel(ODFdiscrete.x));
ODFdiscrete.x = ODFdiscrete.x(RandSample);
ODFdiscrete.y = ODFdiscrete.y(RandSample);
ODFdiscrete.z = ODFdiscrete.z(RandSample);
ODFdiscrete.P = ones(size(ODFdiscrete.x));

%Watson distribution smoothing kernel
ODindex = .05;
kappa = 1/tan(ODindex*pi/2);

[K.Xsmooth,K.Xdiscrete] = ndgrid(ODFsmooth.x,ODFdiscrete.x);
[K.Ysmooth,K.Ydiscrete] = ndgrid(ODFsmooth.y,ODFdiscrete.y);
[K.Zsmooth,K.Zdiscrete] = ndgrid(ODFsmooth.z,ODFdiscrete.z);

Kernel = exp(kappa*(K.Xsmooth.*K.Xdiscrete + ...
    K.Ysmooth.*K.Ydiscrete + K.Zsmooth.*K.Zdiscrete).^2);

clear K

%Smooth ODF amplitude

ODFsmooth.P = Kernel*ODFdiscrete.P;
%ODFsmooth.P = ones(size(ODFsmooth.P));
ODFsmooth.P = ODFsmooth.P/sum(ODFsmooth.P)*ODFsmooth.N;

ODFsmooth.vreccum = repmat(ODFsmooth.P,[1 3]).*[sin(ODFsmooth.theta).*cos(ODFsmooth.phi)...
    sin(ODFsmooth.theta).*sin(ODFsmooth.phi) ...
    cos(ODFsmooth.theta)];
ODFsmooth.c = abs([ODFsmooth.x ODFsmooth.y ODFsmooth.z]);


figure(6), clf
axes('position',[0.02 .02 .5*11/8.5 .98],'FontSize',20)
p = patch('Faces',ODFsmooth.tri,'Vertices',ODFsmooth.vreccum);
axis tight, axis square, axis equal
axis(2.5*[-1 1 -1 1 -1 1])
view(30,30)
set(p,'FaceColor','interp','FaceVertexCData',ODFsmooth.c,...
'EdgeColor','k','LineWidth',1)
set(gca,'LineWidth',5,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[])
%axis off
eval(['print ' DataDir '/' ExpNam '/' num2str(expno) '/ODF -djpeg -r300'])

return
%%
subplot(2,2,2)
plot(theta(1,:),Ptheta,'-')

swh = 2.5e3;
si = 1024;
R2 = 10;
DeltaQ = 550;

nu = swh*linspace(0,1,si)'; nu = nu - nu(si/2+1);

[K.nu, K.theta] = ndgrid(nu, theta(1,:));

K.offset = 2*DeltaQ*(3*cos(K.theta).^2-1)/2;

Kernel2H = R2.^2./(R2.^2 + 4*(K.offset/2 - K.nu).^2)...
    + R2.^2./(R2.^2 + 4*(-K.offset/2 - K.nu).^2);

clear K

I2H = Kernel2H*Ptheta';

figure(6)
subplot(2,2,4)
plot(nu,I2H,'-')
axis tight off

return
indx = histdat.v1.z < 0;
histdat.v1.x(indx) = -histdat.v1.x(indx);
histdat.v1.y(indx) = -histdat.v1.y(indx);
histdat.v1.z(indx) = -histdat.v1.z(indx);

[X,Y] = fSchmidt(histdat.v1.x(mask.S0R2CL),histdat.v1.y(mask.S0R2CL),histdat.v1.z(mask.S0R2CL));
latitude.theta = pi/180*[30:30:150 179];
latitude.phi = linspace(0,2*pi,100);
[latitude.phi,latitude.theta] = ndgrid(latitude.phi,latitude.theta);
latitude.z = cos(latitude.theta);
latitude.x = sin(latitude.theta).*cos(latitude.phi);
latitude.y = sin(latitude.theta).*sin(latitude.phi);
[latitude.X,latitude.Y] = fSchmidt(latitude.x,latitude.y,latitude.z);
longitude.theta = pi/180*linspace(30,180,100);
longitude.phi = pi/180*[30:30:360];
[longitude.theta,longitude.phi] = ndgrid(longitude.theta,longitude.phi);
longitude.z = cos(longitude.theta);
longitude.x = sin(longitude.theta).*cos(longitude.phi);
longitude.y = sin(longitude.theta).*sin(longitude.phi);
[longitude.X,longitude.Y] = fSchmidt(longitude.x,longitude.y,longitude.z);

figure(3), clf
plot(latitude.X,latitude.Y,'k-')
hold on
plot(longitude.X,longitude.Y,'k-')
plot(X,Y,'b.','MarkerSize',0.25)
axis tight equal off
%return

maxz = max(r.z)*1e3;
maxx = max(r.x)*1e3;

nx = length(r.x)/2+1;
ny = length(r.y)/2+1;
nz = length(r.z)/2+1;
%%
%V = Images.R2; isovalue = 1/15e-3;;
%V = 1./Images.R2; V(Images.R2 == 0) = 0; isovalue = 15e-3;
V = zeros(Nz,Ny,Nx); V(mask.S0R2CL) = 1; isovalue = 0.5;
GreyVal = 1; RGBval = [0 0 1];

[Z,Y,X] = ndgrid(r.z,r.y,r.x);

figure(7), clf
hist(reshape(V,numel(V),1),20)
indxY = (length(r.y)/2+1):length(r.y);
indxY = 1:length(r.y);
indxZ = round(.25*length(r.z)):round(.75*length(r.z));
V = V(indxZ,indxY,:);
X = X(indxZ,indxY,:);
Y = Y(indxZ,indxY,:);
Z = Z(indxZ,indxY,:);



fv = isosurface(X,Y,Z,V,isovalue);
fvc = isocaps(X,Y,Z,V,isovalue);

figure(6), clf
axh3 = axes('position',[0 0 1 1]);
view(-20,70)
hold on
p = patch(fv);
p2 = patch(fvc);
set(p,'FaceColor',GreyVal*RGBval,'EdgeColor','none');
set(p2,'FaceColor',GreyVal*RGBval,'EdgeColor','none');
daspect([1 1 1])
axis tight, axis equal
%axis(.51*[-1 1 -1 1 -1 1])
axis off
camlight(20, 5) 
lighting phong
hold off

xlim = get(axh3,'XLim'); ylim = get(axh3,'YLim'); zlim = get(axh3,'ZLim');
return
for n = 1:length(indxZ)
    indxZtemp = 1:(length(indxZ)-n);
    fv = isosurface(X(indxZtemp,:,:),Y(indxZtemp,:,:),Z(indxZtemp,:,:),V(indxZtemp,:,:),isovalue);
    fvc = isocaps(X(indxZtemp,:,:),Y(indxZtemp,:,:),Z(indxZtemp,:,:),V(indxZtemp,:,:),isovalue);
    
    figure(4), clf
    axh3 = axes('position',[0 0 1 1]);
    view(-20,30)
    hold on
    p = patch(fv);
    p2 = patch(fvc);
    set(p,'FaceColor',GreyVal*RGBval,'EdgeColor','none');
    set(p2,'FaceColor',GreyVal*RGBval,'EdgeColor','none');
    daspect([1 1 1])
    axis([xlim ylim zlim])
    axis off
    camlight(20, 5) 
    lighting phong
    hold off
    pause(.1)
end

return
figure(2), clf
colormap(hot(256));
clim = [0 1];
for nz = 1:length(r.z)
    Iplot = squeeze(I(nz,:,:)/Imax);
    imagesc(r.x*1e3,r.y*1e3,Iplot,clim)
    set(gca,'YDir','normal','YTick',[])
    axis tight, axis off, axis square
    pause(.1)
end
