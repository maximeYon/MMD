clear all

wd = cd;

%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/opt/topspin2/data/DT/nmr';
DataDir = '/Users/daniel/Dropbox';

ExpNam = 'Yeast_tripgse'; expno = [52:58 67 69 73 75 77]; expno = [73];

Itick = [.01 .1 1]; ItickLabel = {'0.01', '0.1', '1'};
Dtick = [1e-12 1e-11 1e-10 1e-9 1e-8]; DtickLabel = {'-12', '-11', '-10', '-9', '-8'};
C_iso = [0 0 1];
C_aniso = [1 0 0];
fs = 18;
lw = 3;
papersize = [25.4 14]; %for PowerPoint

Imin = 2e-2;
   
for nexp = 1:length(expno)
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/NMRacqus'])
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/FitDat'])
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/isoaniso1D2D_BSdat'])
    BSdat_ILT = BSdat;
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/TricompBSdat'])
    BSdat_Tri = BSdat;
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/BSdat'])
    %figure(1), clf, semilogy(FitDat.Xin{2},FitDat.Yin{2},'o'), return

    Yin = BSdat_ILT.Yin_2D;
    Ycalc = BSdat_ILT.Ycalc_2D;

    bmax = max(BSdat_ILT.Xin_iso);

    Nb_iso = sqrt(AcqDat.Nb);
    Nb_aniso = sqrt(AcqDat.Nb);
    indx_iso = 1:Nb_iso;
    indx_aniso = 1 + Nb_iso*(0:(Nb_aniso-1));


    figure(1), clf

    left = .15;
    bottom = .3;
    width = .35;
    height = .6;


    axh_Sb1D = axes('position',[left bottom width height]);
    lh_aniso1D = semilogy(BSdat_ILT.Xin_aniso/1e9,BSdat_ILT.Ycalc_aniso1D,'r-');
    hold on
    mh_aniso1D = semilogy(BSdat_ILT.Xin_aniso/1e9,BSdat_ILT.Yin_aniso1D,'ro');
    lh_iso1D = semilogy(BSdat_ILT.Xin_iso/1e9,BSdat_ILT.Ycalc_iso1D,'b-');
    mh_iso1D = semilogy(BSdat_ILT.Xin_iso/1e9,BSdat_ILT.Yin_iso1D,'bo');
    set([mh_iso1D mh_aniso1D],'LineWidth',.1*lw)
    set([lh_iso1D lh_aniso1D],'LineWidth',lw)
    set([lh_iso1D mh_iso1D],'Color',C_iso)
    set([lh_aniso1D mh_aniso1D],'Color',C_aniso)
    set([mh_aniso1D],'MarkerFaceColor',C_aniso,'MarkerSize',5*lw)
    set([mh_iso1D],'MarkerFaceColor',C_iso,'MarkerSize',4*lw)
    set(gca,'YLim',[.04 1.1],'XLim',max(BSdat_ILT.Xin_iso/1e9)*[-.1 1.1])
    set(gca,'YTick',Itick,'YTickLabel',ItickLabel,'TickLength',.01/width*[1 1])
    
    left = 1-width-.04;

    axh_PD1D = axes('position',[left bottom width height]);
    lh_iso1D = semilogx(BSdat_ILT.Diso,BSdat_ILT.PD_iso,'r-');
    hold on
    lh_aniso1D = semilogx(BSdat_ILT.Daniso,BSdat_ILT.PD_aniso + 1.1*max(BSdat_ILT.PD_iso),'b-');
    axis tight
    ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.1 1.05]; set(gca,'YLim',ylim)
    %xlim = get(gca,'XLim'); xlim = [xlim(1) 2*xlim(2)]; set(gca,'XLim',xlim)
    set([lh_iso1D lh_aniso1D],'LineWidth',lw)
    set([lh_iso1D],'Color',C_iso)
    set([lh_aniso1D],'Color',C_aniso)
    set(gca,'TickDir','out','Box','off','LineWidth',lw,'FontSize',fs,'XTick',Dtick,'TickLength',.01/width*[1 1],'YTick',[],'YAxisLocation','right')

    set([axh_Sb1D; axh_PD1D],'TickDir','out','Box','off','LineWidth',lw,'FontSize',fs)

    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    eval(['print -dpdf -loose ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/BackaskogFig1D'])

    figure(2), clf

    bottom = .3;
    left = .15;
    width = .3;
    height = .6;
    
    Xval = 1e-9*reshape(AcqDat.bmat.iso(1:AcqDat.Nb),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    Yval = 1e-9*reshape(AcqDat.bmat.aniso(1:AcqDat.Nb),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));

    axh1 = axes('position',[left bottom width height]);
    Zval = reshape(log10(Ycalc),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    Zval(Zval<log10(Imin)) = NaN;
    
    lh1 = plot3(Xval,Yval,Zval,'k-');
    hold on
    lh2 = plot3(Xval',Yval',Zval','k-');
    lh3 = plot3(Xval(indx_iso),Yval(indx_iso),Zval(indx_iso),'k-');
    lh4 = plot3(Xval(indx_aniso)',Yval(indx_aniso)',Zval(indx_aniso)','k-');

    Zval = reshape(log10(Yin),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    Zval(Zval<log10(Imin)) = NaN;
    mh1 = plot3(Xval,Yval,Zval,'ko');
    mh2 = plot3(Xval(indx_iso),Yval(indx_iso),Zval(indx_iso),'ko');
    mh3 = plot3(Xval(indx_aniso)',Yval(indx_aniso)',Zval(indx_aniso)','ko');
    view(150,15)
    axis([bmax*1e-9*[-.1 1.1 -.1 1.1] log10(Imin) .1])
    set([lh1; lh2],'LineWidth',0.5*lw)
    set([lh3],'LineWidth',2*lw,'Color',C_iso)
    set([lh4],'LineWidth',2*lw,'Color',C_aniso)
    set([mh1],'LineWidth',.1*lw,'MarkerSize',3*lw,'MarkerFaceColor',[0 0 0])
    set([mh2],'LineWidth',.1*lw,'MarkerSize',4*lw,'Color',C_iso,'MarkerFaceColor',C_iso)
    set([mh3],'LineWidth',.1*lw,'MarkerSize',5*lw,'Color',C_aniso,'MarkerFaceColor',C_aniso)
    set(axh1,'ZTick',log10(Itick),'ZTickLabel',ItickLabel)
    set([axh1 ],'LineWidth',lw,'FontSize',fs,'TickDir','out','Box','off','TickLength',.02/width*[1 1])
    


    minclevel = .05;
    maxclevel = .9;
    Nclevel = 5;
    clevels_norm = linspace(minclevel,maxclevel,Nclevel);
    clevels_norm = logspace(log10(minclevel),log10(maxclevel),Nclevel);

    bottom = .35;
    left = .08;
    width = .2;
    height = width*papersize(1)/papersize(2);
    width_proj = .08;
    height_proj = width_proj*papersize(1)/papersize(2);    
    
    left = 1 - width - .15;
    Xval = log10(BSdat.Di);
    Yval = log10(BSdat.Da);
    Zval = reshape(BSdat.PDiDa,numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1);
    ZprojY = sum(Zval,2);
    DiagXY = [min(Xval) max(Xval)];
    
    axh_sparseILT2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    plot(DiagXY,DiagXY,'--k','LineWidth',.5*lw)
    
    axes('position',[left bottom+height+.01 width height_proj]);
    plot(Xval,ZprojX,'k-','LineWidth',lw,'Color',C_iso)
    axis tight off
    ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)
    axes('position',[left-width_proj-.01 bottom width_proj height]);
    plot(ZprojY,Yval,'k-','LineWidth',lw,'Color',C_aniso)
    set(gca,'XDir','reverse')
    axis tight off
    xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)   
    
    set([axh_sparseILT2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out','XTick',log10(Dtick),'XTickLabel',DtickLabel,'YTick',log10(Dtick),...
        'YTickLabel',DtickLabel,'TickLength',.007/width*[1 1])

    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    eval(['print -dpdf -loose ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/BackaskogFig2D'])
%%
    figure(3), clf
        

    left = .15;
    Xval = log10(BSdat.Dpar);
    Yval = log10(BSdat.Dperp);
    Zval = reshape(BSdat.PDparDperp,numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1);
    ZprojY = sum(Zval,2);
    DiagXY = [min(Xval) max(Xval)];

    axh_sparse2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    hold on
    plot(DiagXY,DiagXY,'--k','LineWidth',.5*lw)

    ConfLev = .95;
    Ncomp = 3;
    for ncomp = 1:Ncomp
        eval(['histdat = log10(BSdat_Tri.Dpar_' num2str(ncomp) ');'])
        histdat = sort(histdat);
        Xmin = histdat(floor((1-ConfLev)/2*BSdat_Tri.NBS));
        Xmax = histdat(ceil((1-(1-ConfLev)/2)*BSdat_Tri.NBS));
        Xmean = mean(histdat);
        eval(['histdat = log10(BSdat_Tri.Dperp_' num2str(ncomp) ');'])
        histdat = sort(histdat);
        Ymin = histdat(floor((1-ConfLev)/2*BSdat_Tri.NBS));
        Ymax = histdat(ceil((1-(1-ConfLev)/2)*BSdat_Tri.NBS));
        Ymean = mean(histdat);
%          plot([Xmin; Xmax],[Ymean; Ymean],'-g','LineWidth',1*lw)
%          plot([Xmean; Xmean],[Ymin; Ymax],'-g','LineWidth',1*lw)
         plot(Xmean,Ymean,'x','LineWidth',1*lw,'MarkerSize',5*lw,'Color',[0 .7 0])
    end
    
    axes('position',[left bottom+height+.01 width height_proj]);
    plot(Xval,ZprojX,'k-','LineWidth',lw)
    axis tight off
    ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)
    axes('position',[left-width_proj-.01 bottom width_proj height]);
    plot(ZprojY,Yval,'k-','LineWidth',lw)
    set(gca,'XDir','reverse')
    axis tight off
    xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)   
    
    left = 1 - width - .15;
    Xval = log10(BSdat.Diso);
    Yval = log10(BSdat.Dratio);
    Zval = reshape(BSdat.PDisoDratio,numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1);
    ZprojY = sum(Zval,2);
    axh_ratio2D = axes('position',[left bottom width height]);
    clevels = max(max(Zval))*clevels_norm;
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    set(gca,'LineWidth',lw,'FontSize',fs,'YDir','normal','YAxisLocation','right','TickDir','out','TickLength',.02*[1 1])
    hold on
    plot(DiagXY,[0 0],'--k','LineWidth',.5*lw)
    
    for ncomp = 1:Ncomp
        eval(['histdat = log10(BSdat_Tri.Diso_' num2str(ncomp) ');'])
        histdat = sort(histdat);
        Xmin = histdat(floor((1-ConfLev)/2*BSdat_Tri.NBS));
        Xmax = histdat(ceil((1-(1-ConfLev)/2)*BSdat_Tri.NBS));
        Xmean = mean(histdat);
        eval(['histdat = log10(BSdat_Tri.Dpar_' num2str(ncomp) './BSdat_Tri.Dperp_' num2str(ncomp) ');'])
        histdat = sort(histdat);
        Ymin = histdat(floor((1-ConfLev)/2*BSdat_Tri.NBS));
        Ymax = histdat(ceil((1-(1-ConfLev)/2)*BSdat_Tri.NBS));
        Ymean = mean(histdat);
%          plot([Xmin; Xmax],[Ymean; Ymean],'-g','LineWidth',1*lw)
%          plot([Xmean; Xmean],[Ymin; Ymax],'-g','LineWidth',1*lw)
         plot(Xmean,Ymean,'x','LineWidth',1*lw,'MarkerSize',5*lw,'Color',[0 .7 0])
    end
    
    axes('position',[left bottom+height+.01 width height_proj]);
    plot(Xval,ZprojX,'k-','LineWidth',lw)
    axis tight off
    ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)
    axes('position',[left-width_proj-.01 bottom width_proj height]);
    plot(ZprojY,Yval,'k-','LineWidth',lw)
    set(gca,'XDir','reverse')
    axis tight off
    xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)   
    set([axh_ratio2D],'XTick',log10(Dtick),'XTickLabel',DtickLabel,'YTick',[-3:3])


    %     bottom = .64;
%     Xval = log10(BSdat_ILT.Diso);
%     Yval = log10(BSdat_ILT.Daniso);
%     Zval = reshape(BSdat_ILT.PD_2D,numel(Xval),numel(Yval))';
%     ZprojX = sum(Zval,1);
%     ZprojY = sum(Zval,2);
%     axh_ILT2D = axes('position',[left bottom width width*11/8.5]);
%     clevels = max(max(Zval))*clevels_norm;
%     contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
%     hold on
%     plot(DiagXY,DiagXY,'--k','LineWidth',.5*lw)
%     
%     axes('position',[left bottom+width*11/8.5+.01 width height_proj]);
%     plot(Xval,ZprojX,'k-','LineWidth',lw,'Color',C_iso)
%     axis tight off
%     ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)
%     axes('position',[left-height_proj*8.5/11-.01 bottom height_proj*8.5/11 width*11/8.5]);
%     plot(ZprojY,Yval,'k-','LineWidth',lw,'Color',C_aniso)
%     set(gca,'XDir','reverse')
%     axis tight off
%     xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)   

    set([axh_sparse2D; axh_ratio2D],'LineWidth',lw,'FontSize',fs,...
        'YDir','normal','YAxisLocation','right','TickDir','out',...
        'TickLength',.007/width*[1 1])
    set([axh_sparse2D],'XTick',log10(Dtick),'XTickLabel',DtickLabel,'YTick',log10(Dtick),...
        'YTickLabel',DtickLabel)

    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    eval(['print -dpdf -loose ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/BackaskogFig2D_DTD'])
end

cd(wd)