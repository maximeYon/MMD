clear all

fs = 20;

%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/Users/daniel/Dropbox';
DataDir = '/Users/daniel/Documents/Spaces/Presentations';

ExpNam = 'C14E5'; expno = 324:341;
ExpNam = 'C14E5_2'; expno = [4 7:8];

ExpNam = 'RandSampDR1R2'; expno = 35;


for nexp = 1:length(expno)
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/Images'])
    ImagesDT = Images;
    Nvoxels = numel(Images.S0);
    Nx = length(r.x);
    Ny = length(r.y);
    Nz = length(r.z);

    %%
    histdat.DTS0 = reshape(ImagesDT.S0,Nvoxels,1);
    histdat.ADC = reshape(ImagesDT.ADC,Nvoxels,1);
    histdat.FA = reshape(ImagesDT.FA,Nvoxels,1);
    histdat.CL = reshape(ImagesDT.CL,Nvoxels,1);
    histdat.CP = reshape(ImagesDT.CP,Nvoxels,1);
    histdat.v3.x = reshape(ImagesDT.v3.x,Nvoxels,1);
    histdat.v3.y = reshape(ImagesDT.v3.y,Nvoxels,1);
    histdat.v3.z = reshape(ImagesDT.v3.z,Nvoxels,1);
    histdat.P2 = (3*histdat.v3.z.^2 - 1)/2;
    histdat.theta = acos(sqrt((2*histdat.P2 + 1)/3));
    histdat.z = reshape(repmat(r.z, [1 Ny Nx]),Nvoxels,1);

%     mask = find(histdat.DTS0>0 ...
%         & histdat.CP>0*histdat.CL ...
%         & histdat.z > -6e-3 & histdat.z < 6e-3 ...
%         & histdat.theta > 0/180*pi & histdat.theta < 90/180*pi);
    mask = find(histdat.DTS0>0 ...
        & histdat.CL>1*histdat.CP ...
        & histdat.z > -6e-3 & histdat.z < 6e-3 ...
        & histdat.theta > 0/180*pi & histdat.theta < 90/180*pi);

    Nsubx = 4;
    Nsuby = 2;
    Nbins = 500;
    figure(1), clf
    subplot(Nsuby,Nsubx,2)
    hist(histdat.ADC(mask)*1e9,Nbins)
    xlabel('ADC')
    subplot(Nsuby,Nsubx,3)
    hist(histdat.FA(mask),Nbins)
    xlabel('FA')
    subplot(Nsuby,Nsubx,4)
    hist(histdat.CL(mask),Nbins)
    xlabel('CL')
    subplot(Nsuby,Nsubx,5)
    hist(histdat.CP(mask),Nbins)
    xlabel('CP')
    subplot(Nsuby,Nsubx,6)
    hist(histdat.theta(mask)/pi*180,Nbins)
    xlabel('\theta')


    %%

    Nslices = 8^2;
    slices = round(Nz*linspace(.2,.8,Nslices+2)); slices([1 Nslices+2]) = [];
    width = .5/sqrt(Nslices);
    height = 11/8.5*.5/sqrt(Nslices);

    figure(4), clf           

    c.bright = zeros(Nz,Ny,Nx);
%     c.bright(mask) = ImagesDT.CP(mask);
%     c.r = abs(ImagesDT.v3.x);
%     c.g = abs(ImagesDT.v3.y);
%     c.b = abs(ImagesDT.v3.z);
    c.bright(mask) = ImagesDT.CL(mask);
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

    %eval(['print ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/CP -depsc -loose'])
    eval(['print ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/CP -depsc -loose'])
end

