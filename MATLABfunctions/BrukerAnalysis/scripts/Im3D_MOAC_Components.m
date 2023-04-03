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
    load(['PDisoR2_a'])

    [Nx,Ny,Nz,Nc] = size(PDisoR2_a);
    Nz = 1;

    PDisoR2 = squeeze(sum(sum(sum(PDisoR2_a,3),2),1));

    minclevel = .1;
    maxclevel = .9;
    Nclevel = 5;
    clevels_norm = linspace(minclevel,maxclevel,Nclevel);
    clevels_norm = logspace(log10(minclevel),log10(maxclevel),Nclevel);
%%
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
    Yval = log10(R2);
    Zval = PDisoR2;
    XLabel = [];
    YLabel = [];
    xtick = [-11:.5:8];
    ytick = [0:.5:2];
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
    
    nx = 25; ny = 20; nz = 1;
    nx = 26; ny = 37; nz = 1;
    nx = 55; ny = 48; nz = 1;
    nx = 32; ny = 30; nz = 1;
    nx = 40; ny = 37; nz = 1;
    nx = 36; ny = 11; nz = 1;
    
    Zval = squeeze(sum(sum(sum(PDisoR2_a(nx,ny,nz,:),3),2),1));

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

    
    [Xval_a,Yval_a] = ndgrid(Xval,Yval);
    Xval_v = reshape(Xval_a,[numel(Xval_a) 1]);
    Yval_v = reshape(Yval_a,[numel(Yval_a) 1]);

    indx = 1:numel(Xval_v);
    C_S0 = squeeze(sum(PDisoR2_a(:,:,:,indx),4))'/2000;
    indx = [Xval_v>-8.75 & Yval_v<0.5];
    C_CSF = squeeze(sum(PDisoR2_a(:,:,:,indx),4))'/2000;
    indx = [Xval_v<-8.75 & Yval_v>1];
    C_tissue = squeeze(sum(PDisoR2_a(:,:,:,indx),4))'/2000;
    indx = [Xval_v>-8.75 & Yval_v>1];
    C_artifact = squeeze(sum(PDisoR2_a(:,:,:,indx),4))'/2000;
    indx = [Xval_v>-9.5 & Xval_v<-9.2 & Yval_v>1 & Yval_v<1.5];
    C_GM = squeeze(sum(PDisoR2_a(:,:,:,indx),4))'/2000;
    indx = [Xval_v>-9.2 & Xval_v<-8.9 & Yval_v>1 & Yval_v<1.5];
    C_WM = squeeze(sum(PDisoR2_a(:,:,:,indx),4))'/2000;
    
    width = .13; height = width*11/8.5;
    bottom = .05; dheight = height+.05;
    left = .05; dleft = width;
    X = (1:Nx)';
    Y = (1:Ny)';

    figure(3), clf
    axes('position',[0 .5 1/3 .5])
    colormap('gray')
    imagesc(X,Y,C_S0)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    hold on
    plot(nx,ny,'ro')
    title('S_0')

    axes('position',[1/3 .5 1/3 .5])
    colormap('gray')
    imagesc(X,Y,C_CSF)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    title('CSF')

    axes('position',[2/3 .5 1/3 .5])
    colormap('gray')
    imagesc(X,Y,C_tissue)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    title('tissue')

    axes('position',[0/3 0 1/3 .5])
    colormap('gray')
    imagesc(X,Y,C_artifact)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    title('artifact')

    axes('position',[1/3 0 1/3 .5])
    colormap('gray')
    imagesc(X,Y,C_GM)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    title('GM')

    axes('position',[2/3 0 1/3 .5])
    colormap('gray')
    imagesc(X,Y,C_WM)
    set(gca,'YDir','normal','Clim',[0 1])
    axis equal, axis tight, axis off
    title('WM')

    set(gcf, 'PaperPosition', [0 0 11 8.5],'PaperSize', [11 8.5]); 
    eval(['print Segmented -loose -dpdf']) 
    
%    save('PDisoR2_a','PDisoR2_a','Diso','R2')    
end
