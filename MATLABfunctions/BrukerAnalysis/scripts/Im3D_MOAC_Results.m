clear all

wd = cd;
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/Users/daniel/Documents/Spaces/Presentations';
%DataDir = '/Users/daniel/Dropbox';
DataDir = '/Users/daniel/Desktop';

ExpNam = {'COMPILED_DATA'};

fs = 12;
lw = 2;

PlotInterm = 0;

cd(DataDir)
for ndir = 1:length(ExpNam)
    cd(ExpNam{ndir})
    load(['ImagesRaw'])
    load(['ImagesMOAC_a'])

    [Nx,Ny,Nz,td1] = size(ImagesRaw.S);
    bT = ImagesRaw.bT;
    Nz = 1;
    ImagesRaw.S = ImagesRaw.S(:,:,15,:);
    Smax = max(reshape(ImagesRaw.S,[numel(ImagesRaw.S) 1]));
    figure(1), clf
    bottom = .7;
    left = .05;
    width = .1;
    height = .1;
    axh_hist1 = axes('position',[left bottom width height]);
    histdat = sqrt(ImagesMOAC_a.chisq); histdat = reshape(histdat,[numel(histdat),1]);
    hist(histdat)
    axis tight
    hold on
    ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
    xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
    plot([xlim'],[0; 0],'k-')
    set(gca,'XLim',xlim,'YLim',ylim)
    xlabel('rms\chi')

    MDmax = 5e-9;

    ImagesMOAC_BS.S0 = sum(ImagesMOAC_a.w,4);
    ImagesMOAC_BS.MD = sum(ImagesMOAC_a.w.*ImagesMOAC_a.iso,4)./(ImagesMOAC_BS.S0+eps);
    ImagesMOAC_BS.VMD = sum(ImagesMOAC_a.w...
        .*(ImagesMOAC_a.iso - repmat(ImagesMOAC_BS.MD,[1 1 1 ReconDat.Nnodes 1])).^2,4)...
    ./(ImagesMOAC_BS.S0+eps);
    ImagesMOAC_BS.R2 = sum(ImagesMOAC_a.w.*ImagesMOAC_a.R2,4)./(ImagesMOAC_BS.S0+eps);
    ImagesMOAC_BS.VR2 = sum(ImagesMOAC_a.w...
        .*(ImagesMOAC_a.R2 - repmat(ImagesMOAC_BS.R2,[1 1 1 ReconDat.Nnodes 1])).^2,4)...
    ./(ImagesMOAC_BS.S0+eps);

    ImagesMOAC.S0 = squeeze(mean(ImagesMOAC_BS.S0,5));
    ImagesMOAC.MD = squeeze(mean(ImagesMOAC_BS.MD,5));
    ImagesMOAC.VMD = squeeze(mean(ImagesMOAC_BS.VMD,5));
    ImagesMOAC.R2 = squeeze(mean(ImagesMOAC_BS.R2,5));
    ImagesMOAC.VR2 = squeeze(mean(ImagesMOAC_BS.VR2,5));

    width = .13; height = width*11/8.5;
    bottom = .05; dheight = height+.05;
    left = .05; dleft = width;
    X = (1:Nx)';
    Y = (1:Ny)';
    C = zeros(numel(Y),numel(X),3);
%%
    figure(3), clf
    axes('position',[0 .5 1/3 .5])
    colormap('gray')
    C = ImagesMOAC.S0'/3000;
    imagesc(X,Y,C)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    title('S_0')
    axes('position',[1/3 .5 1/3 .5])
    colormap('gray')
    C = ImagesMOAC.MD'/3e-9;
    imagesc(X,Y,C)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    title('MD')
    axes('position',[1/3 0 1/3 .5])
    colormap('gray')
    C = (ImagesMOAC.VMD./ImagesMOAC.MD.^2)';
    imagesc(X,Y,C)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    title('VMD/MD^2')
    axes('position',[2/3 .5 1/3 .5])
    colormap('gray')
    C = ImagesMOAC.R2'/20;
    imagesc(X,Y,C)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    title('R2')
    axes('position',[2/3 0 1/3 .5])
    colormap('gray')
    C = 2*(ImagesMOAC.VR2./ImagesMOAC.R2.^2)';
    imagesc(X,Y,C)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    title('VR2/R2^2')
    
    pause(1)

    bottom1 = 0;
    left1 = 0;
    width = 1/Ny;
    height = 1/Nx;

    ND = 32;
    logstd = .02;
    Dmin = ReconDat.Dmin/10^(10*logstd);
    Dmax = ReconDat.Dmax*10^(10*logstd);
    R1min = ReconDat.R1min/10^(10*logstd);
    R1max = ReconDat.R1max*10^(10*logstd);
    R2min = ReconDat.R2min/10^(10*logstd);
    R2max = ReconDat.R2max*10^(10*logstd);

    Disotick = [5e-10 1e-9 2e-9 5e-9];
    Dratiotick = [1e-2 1e-1 1 10 100];
    R1tick = [1e-1 1 10];
    R2tick = [2 5 10 20];

    minclevel = .1;
    maxclevel = .9;
    Nclevel = 5;
    clevels_norm = linspace(minclevel,maxclevel,Nclevel);
    clevels_norm = logspace(log10(minclevel),log10(maxclevel),Nclevel);

    Diso = logspace(log10(Dmin),log10(Dmax),ND);
    Dratio = logspace(log10(Dmin/Dmax),log10(Dmax/Dmin),ND);
    R1 = logspace(log10(R1min),log10(R1max),ND);
    R2 = logspace(log10(R2min),log10(R2max),ND);        

    PDisoDiso_a = zeros(Nx,Ny,Nz,ND^2);
    PDratioDratio_a = zeros(Nx,Ny,Nz,ND^2);
    PR1R1_a = zeros(Nx,Ny,Nz,ND^2);
    PR2R2_a = zeros(Nx,Ny,Nz,ND^2);
    PDisoDratio_a = zeros(Nx,Ny,Nz,ND^2);
    PDisoR1_a = zeros(Nx,Ny,Nz,ND^2);
    PDisoR2_a = zeros(Nx,Ny,Nz,ND^2);
    PDratioR1_a = zeros(Nx,Ny,Nz,ND^2);
    PDratioR2_a = zeros(Nx,Ny,Nz,ND^2);
    PR1R2_a = zeros(Nx,Ny,Nz,ND^2);


    options = optimset('MaxFunEvals',1e4,'Display','off');
    p =  TimedProgressBar( Ny, 10, ...
    'Computing. Remaining time: ', ', Completed: ', 'Concluded in ' );


    for nz = 1:Nz    
        for ny = 1:Ny    
            for nx = 1:Nx
%                     nx = 55; ny = 45; nz = 1;
%                     %nx = 50; ny = 57; nz = 1;
%                     nx = 59; ny = 47; nz = 1;

                uDTrec_a.w = squeeze(ImagesMOAC_a.w(nx,ny,nz,:,:));
                if sum(uDTrec_a.w) > 0
                    uDTrec_a.iso = squeeze(ImagesMOAC_a.iso(nx,ny,nz,:,:));
                    uDTrec_a.Delta = squeeze(ImagesMOAC_a.Delta(nx,ny,nz,:,:));
                    uDTrec_a.theta = squeeze(ImagesMOAC_a.theta(nx,ny,nz,:,:));
                    uDTrec_a.phi = squeeze(ImagesMOAC_a.phi(nx,ny,nz,:,:));
                    uDTrec_a.R1 = squeeze(ImagesMOAC_a.R1(nx,ny,nz,:,:));
                    uDTrec_a.R2 = squeeze(ImagesMOAC_a.R2(nx,ny,nz,:,:));

                    uDTrec.N = ReconDat.NBS;
                    mask = find(uDTrec_a.w>0);
                    uDTrec.w = reshape(uDTrec_a.w(mask),[numel(mask) 1])/uDTrec.N;
                    uDTrec.iso = reshape(uDTrec_a.iso(mask),[numel(mask) 1]);
                    uDTrec.Delta = reshape(uDTrec_a.Delta(mask),[numel(mask) 1]);
                    uDTrec.theta = reshape(uDTrec_a.theta(mask),[numel(mask) 1]);
                    uDTrec.phi = reshape(uDTrec_a.phi(mask),[numel(mask) 1]);
                    uDTrec.R1 = reshape(uDTrec_a.R1(mask),[numel(mask) 1]);
                    uDTrec.R2 = reshape(uDTrec_a.R2(mask),[numel(mask) 1]);
                    uDTrec.par = uDTrec.iso.*(1 + 2*uDTrec.Delta);
                    uDTrec.perp = uDTrec.iso.*(1 - uDTrec.Delta);
                    uDTrec.ratio = (1 + 2*uDTrec.Delta)./(1 - uDTrec.Delta);
                    uDTrec.N = numel(mask);

                    [uDTrec.w,indx] = sort(uDTrec.w,'descend');
                    uDTrec.iso = uDTrec.iso(indx);
                    uDTrec.Delta = uDTrec.Delta(indx);
                    uDTrec.theta = uDTrec.theta(indx);
                    uDTrec.phi = uDTrec.phi(indx);
                    uDTrec.R1 = uDTrec.R1(indx);
                    uDTrec.R2 = uDTrec.R2(indx);
                    uDTrec.par = uDTrec.par(indx);
                    uDTrec.perp = uDTrec.iso(indx);
                    uDTrec.ratio = uDTrec.ratio(indx);


                    X = Diso; Y = Diso; Xrec = uDTrec.iso; Yrec = uDTrec.iso;

                    [X_a,Y_a] = ndgrid(X,Y);
                    X_v = reshape(X_a,ND^2,1);
                    Y_v = reshape(Y_a,ND^2,1);
                    [K1,K1u] = ndgrid(X_v,Xrec);
                    [K2,K2u] = ndgrid(Y_v,Yrec);
                    Kernel = exp(-((log10(K1) - log10(K1u)).^2....
                    + (log10(K2) - log10(K2u)).^2)/2/logstd^2);
                    PD = Kernel*uDTrec.w;

                    PDisoDiso = sum(uDTrec.w)*PD/sum(PD);

                    X = Dratio; Y = Dratio; Xrec = uDTrec.ratio; Yrec = uDTrec.ratio;

                    [X_a,Y_a] = ndgrid(X,Y);
                    X_v = reshape(X_a,ND^2,1);
                    Y_v = reshape(Y_a,ND^2,1);
                    [K1,K1u] = ndgrid(X_v,Xrec);
                    [K2,K2u] = ndgrid(Y_v,Yrec);
                    Kernel = exp(-((log10(K1) - log10(K1u)).^2)/2/(2*logstd)^2).*...
                        exp(-((log10(K2) - log10(K2u)).^2)/2/(2*logstd)^2);
                    PD = Kernel*uDTrec.w;

                    PDratioDratio = sum(uDTrec.w)*PD/sum(PD);

                    X = R1; Y = R1; Xrec = uDTrec.R1; Yrec = uDTrec.R1;

                    [X_a,Y_a] = ndgrid(X,Y);
                    X_v = reshape(X_a,ND^2,1);
                    Y_v = reshape(Y_a,ND^2,1);
                    [K1,K1u] = ndgrid(X_v,Xrec);
                    [K2,K2u] = ndgrid(Y_v,Yrec);
                    Kernel = exp(-((log10(K1) - log10(K1u)).^2....
                    + (log10(K2) - log10(K2u)).^2)/2/logstd^2);
                    PD = Kernel*uDTrec.w;

                    PR1R1 = sum(uDTrec.w)*PD/sum(PD);

                    X = R2; Y = R2; Xrec = uDTrec.R2; Yrec = uDTrec.R2;

                    [X_a,Y_a] = ndgrid(X,Y);
                    X_v = reshape(X_a,ND^2,1);
                    Y_v = reshape(Y_a,ND^2,1);
                    [K1,K1u] = ndgrid(X_v,Xrec);
                    [K2,K2u] = ndgrid(Y_v,Yrec);
                    Kernel = exp(-((log10(K1) - log10(K1u)).^2....
                    + (log10(K2) - log10(K2u)).^2)/2/logstd^2);
                    PD = Kernel*uDTrec.w;

                    PR2R2 = sum(uDTrec.w)*PD/sum(PD);

                    X = Diso; Y = Dratio; Xrec = uDTrec.iso; Yrec = uDTrec.ratio;

                    [X_a,Y_a] = ndgrid(X,Y);
                    X_v = reshape(X_a,ND^2,1);
                    Y_v = reshape(Y_a,ND^2,1);
                    [K1,K1u] = ndgrid(X_v,Xrec);
                    [K2,K2u] = ndgrid(Y_v,Yrec);
                    Kernel = exp(-((log10(K1) - log10(K1u)).^2)/2/(1*logstd)^2).*...
                        exp(-((log10(K2) - log10(K2u)).^2)/2/(2*logstd)^2);
                    PD = Kernel*uDTrec.w;

                    PDisoDratio = sum(uDTrec.w)*PD/sum(PD);

                    X = Diso; Y = R1; Xrec = uDTrec.iso; Yrec = uDTrec.R1;

                    [X_a,Y_a] = ndgrid(X,Y);
                    X_v = reshape(X_a,ND^2,1);
                    Y_v = reshape(Y_a,ND^2,1);
                    [K1,K1u] = ndgrid(X_v,Xrec);
                    [K2,K2u] = ndgrid(Y_v,Yrec);
                    Kernel = exp(-((log10(K1) - log10(K1u)).^2....
                    + (log10(K2) - log10(K2u)).^2)/2/logstd^2);
                    PD = Kernel*uDTrec.w;

                    PDisoR1 = sum(uDTrec.w)*PD/sum(PD);

                    X = Diso; Y = R2; Xrec = uDTrec.iso; Yrec = uDTrec.R2;

                    [X_a,Y_a] = ndgrid(X,Y);
                    X_v = reshape(X_a,ND^2,1);
                    Y_v = reshape(Y_a,ND^2,1);
                    [K1,K1u] = ndgrid(X_v,Xrec);
                    [K2,K2u] = ndgrid(Y_v,Yrec);
                    Kernel = exp(-((log10(K1) - log10(K1u)).^2....
                    + (log10(K2) - log10(K2u)).^2)/2/logstd^2);
                    PD = Kernel*uDTrec.w;

                    PDisoR2 = sum(uDTrec.w)*PD/sum(PD);

                    X = Dratio; Y = R1; Xrec = uDTrec.ratio; Yrec = uDTrec.R1;

                    [X_a,Y_a] = ndgrid(X,Y);
                    X_v = reshape(X_a,ND^2,1);
                    Y_v = reshape(Y_a,ND^2,1);
                    [K1,K1u] = ndgrid(X_v,Xrec);
                    [K2,K2u] = ndgrid(Y_v,Yrec);
                    Kernel = exp(-((log10(K1) - log10(K1u)).^2)/2/(2*logstd)^2).*...
                        exp(-((log10(K2) - log10(K2u)).^2)/2/(1*logstd)^2);
                    PD = Kernel*uDTrec.w;

                    PDratioR1 = sum(uDTrec.w)*PD/sum(PD);

                    X = Dratio; Y = R2; Xrec = uDTrec.ratio; Yrec = uDTrec.R2;

                    [X_a,Y_a] = ndgrid(X,Y);
                    X_v = reshape(X_a,ND^2,1);
                    Y_v = reshape(Y_a,ND^2,1);
                    [K1,K1u] = ndgrid(X_v,Xrec);
                    [K2,K2u] = ndgrid(Y_v,Yrec);
                    Kernel = exp(-((log10(K1) - log10(K1u)).^2)/2/(2*logstd)^2).*...
                        exp(-((log10(K2) - log10(K2u)).^2)/2/(1*logstd)^2);
                    PD = Kernel*uDTrec.w;

                    PDratioR2 = sum(uDTrec.w)*PD/sum(PD);

                    X = R1; Y = R2; Xrec = uDTrec.R1; Yrec = uDTrec.R2;

                    [X_a,Y_a] = ndgrid(X,Y);
                    X_v = reshape(X_a,ND^2,1);
                    Y_v = reshape(Y_a,ND^2,1);
                    [K1,K1u] = ndgrid(X_v,Xrec);
                    [K2,K2u] = ndgrid(Y_v,Yrec);
                    Kernel = exp(-((log10(K1) - log10(K1u)).^2....
                    + (log10(K2) - log10(K2u)).^2)/2/logstd^2);
                    PD = Kernel*uDTrec.w;

                    PR1R2 = sum(uDTrec.w)*PD/sum(PD);

                    width = .15;
                    height = width*11/8.5;
                    height_proj = .08;
                    ticklength = .007/width*[1 1];
                    dbottom = height + .01;
                    dleft = width + .01;
                    bottom1 = .1;
                    left1 = .1;

                    if PlotInterm
                        Xval = log10(Diso);
                        Yval = log10(Diso);
                        Zval = PDisoDiso;
                        XLabel = 'log(D_{iso}/m^2s^-^1)';
                        YLabel = [];
                        xtick = log10(Disotick);
                        ytick = [];
                        left = left1;
                        bottom = bottom1;

                        figure(1), clf

                        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)
                        axes('position',[left-height_proj*8.5/11-.01 bottom height_proj*8.5/11 width*11/8.5]);
                        plot(ZprojY,Yval,'k-','LineWidth',lw)
                        hold on
                        set(gca,'XDir','reverse')
                        axis tight off
                        xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(Diso);
                        Yval = log10(Dratio);
                        Zval = PDisoDratio;
                        XLabel = [];
                        YLabel = [];
                        xtick = [];
                        ytick = [];
                        left = left1;
                        bottom = bottom1  + dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)
                        axes('position',[left-height_proj*8.5/11-.01 bottom height_proj*8.5/11 width*11/8.5]);
                        plot(ZprojY,Yval,'k-','LineWidth',lw)
                        hold on
                        set(gca,'XDir','reverse')
                        axis tight off
                        xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)


                        Xval = log10(Diso);
                        Yval = log10(R1);
                        Zval = PDisoR1;
                        XLabel = [];
                        YLabel = [];
                        xtick = [];
                        ytick = [];
                        left = left1;
                        bottom = bottom1  + 2*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)
                        axes('position',[left-height_proj*8.5/11-.01 bottom height_proj*8.5/11 width*11/8.5]);
                        plot(ZprojY,Yval,'k-','LineWidth',lw)
                        hold on
                        set(gca,'XDir','reverse')
                        axis tight off
                        xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)


                        Xval = log10(Diso);
                        Yval = log10(R2);
                        Zval = PDisoR2;
                        XLabel = [];
                        YLabel = [];
                        xtick = [];
                        ytick = [];
                        left = left1;
                        bottom = bottom1  + 3*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)
                        axes('position',[left bottom+width*11/8.5+.01 width height_proj]);
                        plot(Xval,ZprojX,'k-','LineWidth',lw)
                        hold on
                        axis tight off
                        ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)
                        axes('position',[left-height_proj*8.5/11-.01 bottom height_proj*8.5/11 width*11/8.5]);
                        plot(ZprojY,Yval,'k-','LineWidth',lw)
                        hold on
                        set(gca,'XDir','reverse')
                        axis tight off
                        xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(Dratio);
                        Yval = log10(R2);
                        Zval = PDratioR2;
                        XLabel = [];
                        YLabel = [];
                        xtick = [];
                        ytick = [];
                        left = left1 + dleft;
                        bottom = bottom1  + 3*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)
                        axes('position',[left bottom+width*11/8.5+.01 width height_proj]);
                        plot(Xval,ZprojX,'k-','LineWidth',lw)
                        hold on
                        axis tight off
                        ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(R1);
                        Yval = log10(R2);
                        Zval = PR1R2;
                        XLabel = [];
                        YLabel = [];
                        xtick = [];
                        ytick = [];
                        left = left1 + 2*dleft;
                        bottom = bottom1  + 3*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)
                        axes('position',[left bottom+width*11/8.5+.01 width height_proj]);
                        plot(Xval,ZprojX,'k-','LineWidth',lw)
                        hold on
                        axis tight off
                        ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(R2);
                        Yval = log10(R2);
                        Zval = PR2R2;
                        XLabel = [];
                        YLabel = 'log(R_{2} / s^-^1)';
                        xtick = [];
                        ytick = log10(R2tick);
                        left = left1 + 3*dleft;
                        bottom = bottom1  + 3*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)
                        axes('position',[left bottom+width*11/8.5+.01 width height_proj]);
                        plot(Xval,ZprojX,'k-','LineWidth',lw)
                        hold on
                        axis tight off
                        ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(Dratio);
                        Yval = log10(Diso);
                        Zval = PDisoDratio;
                        XLabel = 'log(D_{||}/D_\perp)';
                        YLabel = [];
                        xtick = log10(Dratiotick);
                        ytick = [];
                        left = left1 + 1*dleft;
                        bottom = bottom1  + 0*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
                %        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(R1);
                        Yval = log10(Diso);
                        Zval = PDisoR1;
                        XLabel = 'log(R_{1} / s^-^1)';
                        YLabel = [];
                        xtick = log10(R1tick);
                        ytick = [];
                        left = left1 + 2*dleft;
                        bottom = bottom1  + 0*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
                %        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)


                        Xval = log10(R1);
                        Yval = log10(Dratio);
                        Zval = PDratioR1;
                        XLabel = [];
                        YLabel = [];
                        xtick = [];
                        ytick = [];
                        left = left1 + 2*dleft;
                        bottom = bottom1  + 1*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
                %        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(R2);
                        Yval = log10(Dratio);
                        Zval = PDratioR2;
                        XLabel = [];
                        YLabel = 'log(D_{||}/D_\perp)';
                        xtick = [];
                        ytick = log10(Dratiotick);
                        left = left1 + 3*dleft;
                        bottom = bottom1  + 1*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
                %        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(R2);
                        Yval = log10(Diso);
                        Zval = PDisoR2;
                        XLabel = 'log(R_{2} / s^-^1)';
                        YLabel = 'log(D_{iso}/m^2s^-^1)';
                        xtick = log10(R2tick);
                        ytick = log10(Disotick);
                        left = left1 + 3*dleft;
                        bottom = bottom1  + 0*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
                %        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(R2);
                        Yval = log10(R1);
                        Zval = PR1R2;
                        XLabel = [];
                        YLabel = 'log(R_{1} / s^-^1)';
                        xtick = [];
                        ytick = log10(R1tick);
                        left = left1 + 3*dleft;
                        bottom = bottom1  + 2*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
                %        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(Dratio);
                        Yval = log10(Dratio);
                        Zval = PDratioDratio;
                        XLabel = [];
                        YLabel = [];
                        xtick = [];
                        ytick = [];
                        left = left1 + 1*dleft;
                        bottom = bottom1  + 1*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(R1);
                        Yval = log10(R1);
                        Zval = PR1R1;
                        XLabel = [];
                        YLabel = [];
                        xtick = [];
                        ytick = [];
                        left = left1 + 2*dleft;
                        bottom = bottom1  + 2*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        Xval = log10(Dratio);
                        Yval = log10(R1);
                        Zval = PDratioR1;
                        XLabel = [];
                        YLabel = [];
                        xtick = [];
                        ytick = [];
                        left = left1 + 1*dleft;
                        bottom = bottom1  + 2*dbottom;

                        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
                        ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
                        ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

                        axh_2D = axes('position',[left bottom width height]);
                        clevels = max(max(Zval))*clevels_norm;
                        contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
                        hold on
                        xlabel(XLabel), ylabel(YLabel)

                        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
                            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
                            'TickLength',ticklength)

                        figure(4), clf
                        axes('position',[0 0 1 1])
                        colormap('gray')
                        X = (1:Nx)';
                        Y = (1:Ny)';
                        C = ImagesMOAC.S0'/3000;
                        imagesc(X,Y,C)
                        hold on
                        plot(nx,ny,'ro')
                        set(gca,'YDir','normal','Clim',[0 1])
                        axis equal, axis tight, axis off
                        title('S_0')
                        pause(.1)
                    end

                    PDisoDiso_a(nx,ny,nz,:) = PDisoDiso;
                    PDratioDratio_a(nx,ny,nz,:) = PDratioDratio;
                    PR1R1_a(nx,ny,nz,:) = PR1R1;
                    PR2R2_a(nx,ny,nz,:) = PR2R2;
                    PDisoDratio_a(nx,ny,nz,:) = PDisoDratio;
                    PDisoR1_a(nx,ny,nz,:) = PDisoR1;
                    PDisoR2_a(nx,ny,nz,:) = PDisoR2;
                    PDratioR1_a(nx,ny,nz,:) = PDratioR1;
                    PDratioR2_a(nx,ny,nz,:) = PDratioR2;
                    PR1R2_a(nx,ny,nz,:) = PR1R2;

                end
            end
            p.progress; %Counter for progress report
        end
        p.stop;
    end
%%
    PDisoDiso = squeeze(sum(sum(sum(PDisoDiso_a,3),2),1));
    PDratioDratio = squeeze(sum(sum(sum(PDratioDratio_a,3),2),1));
    PR1R1 = squeeze(sum(sum(sum(PR1R1_a,3),2),1));
    PR2R2 = squeeze(sum(sum(sum(PR2R2_a,3),2),1));
    PDisoDratio = squeeze(sum(sum(sum(PDisoDratio_a,3),2),1));
    PDisoR1 = squeeze(sum(sum(sum(PDisoR1_a,3),2),1));
    PDisoR2 = squeeze(sum(sum(sum(PDisoR2_a,3),2),1));
    PDratioR1 = squeeze(sum(sum(sum(PDratioR1_a,3),2),1));
    PDratioR2 = squeeze(sum(sum(sum(PDratioR2_a,3),2),1));
    PR1R2 = squeeze(sum(sum(sum(PR1R2_a,3),2),1));
    
    figure(1), clf

    width = .19;
    height = width*1;
    height_proj = .08;
    ticklength = .007/width*[1 1];
    dbottom = height + .01;
    dleft = width + .01;
    bottom1 = .08;
    left1 = .11;

    Xval = log10(Diso);
    Yval = log10(Diso);
    Zval = PDisoDiso;
    XLabel = [];
    YLabel = [];
    xtick = [];
    ytick = [];
    left = left1 + 2*dleft;
    bottom = bottom1 + 2*dbottom;


    Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(Diso);
    Yval = log10(Dratio);
    Zval = PDisoDratio;
    XLabel = [];
    YLabel = [];
    xtick = [];
    ytick = [];
    left = left1 + 2*dleft;
    bottom = bottom1  + 3*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)
    axes('position',[left bottom+width*1+.01 width height_proj]);
    plot(Xval,ZprojX,'k-','LineWidth',lw)
    hold on
    axis tight off
    ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)


    Xval = log10(Diso);
    Yval = log10(R1);
    Zval = PDisoR1;
    XLabel = 'log(D_{iso}/m^2s^-^1)';
    YLabel = [];
    xtick = log10(Disotick);
    ytick = [];
    left = left1 + 2*dleft;
    bottom = bottom1  + 0*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)


    Xval = log10(Diso);
    Yval = log10(R2);
    Zval = PDisoR2;
    XLabel = [];
    YLabel = [];
    xtick = [];
    ytick = [];
    left = left1 + 2*dleft;
    bottom = bottom1  + 1*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(Dratio);
    Yval = log10(R2);
    Zval = PDratioR2;
    XLabel = [];
    YLabel = 'log(R_{2} / s^-^1)';
    xtick = [];
    ytick = log10(R2tick);
    left = left1 + 3*dleft;
    bottom = bottom1  + 1*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(R1);
    Yval = log10(R2);
    Zval = PR1R2;
    XLabel = [];
    YLabel = [];
    xtick = [];
    ytick = [];
    left = left1 + 0*dleft;
    bottom = bottom1  + 1*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)
    axes('position',[left-height_proj*1-.01 bottom height_proj*1 width*1]);
    plot(ZprojY,Yval,'k-','LineWidth',lw)
    hold on
    set(gca,'XDir','reverse')
    axis tight off
    xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(R2);
    Yval = log10(R2);
    Zval = PR2R2;
    XLabel = [];
    YLabel = [];
    xtick = [];
    ytick = [];
    left = left1 + 1*dleft;
    bottom = bottom1  + 1*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(Dratio);
    Yval = log10(Diso);
    Zval = PDisoDratio;
    XLabel = [];
    YLabel = 'log(D_{iso}/m^2s^-^1)';
    xtick = [];
    ytick = log10(Disotick);
    left = left1 + 3*dleft;
    bottom = bottom1  + 2*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
%        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(R1);
    Yval = log10(Diso);
    Zval = PDisoR1;
    XLabel = [];
    YLabel = [];
    xtick = [];
    ytick = [];
    left = left1 + 0*dleft;
    bottom = bottom1  + 2*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
%        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)
    axes('position',[left-height_proj*1-.01 bottom height_proj*1 width*1]);
    plot(ZprojY,Yval,'k-','LineWidth',lw)
    hold on
    set(gca,'XDir','reverse')
    axis tight off
    xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)


    Xval = log10(R1);
    Yval = log10(Dratio);
    Zval = PDratioR1;
    XLabel = [];
    YLabel = [];
    xtick = [];
    ytick = [];
    left = left1 + 0*dleft;
    bottom = bottom1  + 3*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
%        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)
    axes('position',[left-height_proj*1-.01 bottom height_proj*1 width*1]);
    plot(ZprojY,Yval,'k-','LineWidth',lw)
    hold on
    set(gca,'XDir','reverse')
    axis tight off
    xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    
    axes('position',[left bottom+width*1+.01 width height_proj]);
    plot(Xval,ZprojX,'k-','LineWidth',lw)
    hold on
    axis tight off
    ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(R2);
    Yval = log10(Dratio);
    Zval = PDratioR2;
    XLabel = [];
    YLabel = [];
    xtick = [];
    ytick = [];
    left = left1 + 1*dleft;
    bottom = bottom1  + 3*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
%        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)
    axes('position',[left bottom+width*1+.01 width height_proj]);
    plot(Xval,ZprojX,'k-','LineWidth',lw)
    hold on
    axis tight off
    ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(R2);
    Yval = log10(Diso);
    Zval = PDisoR2;
    XLabel = [];
    YLabel = [];
    xtick = [];
    ytick = [];
    left = left1 + 1*dleft;
    bottom = bottom1  + 2*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
%        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(R2);
    Yval = log10(R1);
    Zval = PR1R2;
    XLabel = 'log(R_{2} / s^-^1)';
    YLabel = [];
    xtick = log10(R2tick);
    ytick = [];
    left = left1 + 1*dleft;
    bottom = bottom1  + 0*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Yval),numel(Xval));
%        Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(Dratio);
    Yval = log10(Dratio);
    Zval = PDratioDratio;
    XLabel = [];
    YLabel = 'log(D_{||}/D_\perp)';
    xtick = [];
    ytick = log10(Dratiotick);
    left = left1 + 3*dleft;
    bottom = bottom1  + 3*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)
    axes('position',[left bottom+width*1+.01 width height_proj]);
    plot(Xval,ZprojX,'k-','LineWidth',lw)
    hold on
    axis tight off
    ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(R1);
    Yval = log10(R1);
    Zval = PR1R1;
    XLabel = 'log(R_{1} / s^-^1)';
    YLabel = [];
    xtick = log10(R1tick);
    ytick = [];
    left = left1 + 0*dleft;
    bottom = bottom1  + 0*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)
    axes('position',[left-height_proj*1-.01 bottom height_proj*1 width*1]);
    plot(ZprojY,Yval,'k-','LineWidth',lw)
    hold on
    set(gca,'XDir','reverse')
    axis tight off
    xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    Xval = log10(Dratio);
    Yval = log10(R1);
    Zval = PDratioR1;
    XLabel = 'log(D_{||}/D_\perp)';
    YLabel = 'log(R_{1} / s^-^1)';
    xtick = log10(Dratiotick);
    ytick = log10(R1tick);
    left = left1 + 3*dleft;
    bottom = bottom1  + 0*dbottom;

    Zval = reshape(Zval/max(Zval),numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1); ZprojX = ZprojX/max(ZprojX);
    ZprojY = sum(Zval,2); ZprojY = ZprojY/max(ZprojY);

    axh_2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    xlabel(XLabel), ylabel(YLabel)

    set([axh_2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
        'TickLength',ticklength)

    set(gcf, 'PaperPosition', [0 0 20 20],'PaperSize', [20 20]); 
    eval(['print MOACFig -loose -dpdf']) 
    
    save('PDisoR2_a','PDisoR2_a','Diso','R2')    
end
