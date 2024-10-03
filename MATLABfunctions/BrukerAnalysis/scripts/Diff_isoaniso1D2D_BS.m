clear all

wd = cd;

%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/opt/topspin2/data/DT/nmr';
DataDir = '/Users/daniel/Dropbox';

ExpNam = 'Yeast_tripgse'; expno = [52:58 67 69 73 75 77]; expno = [73 75 77]; expno = [75 77];

td1start = 2;
signal = 'area';
PlotInterm = 0;
NBS = 1e4;
Dmin = 1e-12; Dmax = 10e-9;
NDmin = 50; NDmax = 100;
logstd = .05; %standard deviation for lognormal PD smoothing
ND_smooth = 64;
NDtot_smooth = ND_smooth^2;

for nexp = 1:length(expno)
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/NMRacqus'])
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/FitDat'])
    %figure(1), clf, semilogy(FitDat.Xin{2},FitDat.Yin{2},'o'), return

    waterpeak = PP.Npeaks;
    PP.peakppm(PP.td1start,waterpeak)
    S = PP.Apeak(:,waterpeak);
    if isfield(AcqDat,'NDeltab')
        AcqDat.Nb = AcqDat.Nb*AcqDat.NDeltab;
    end
    if isfield(AcqDat,'Nvd')==0
        AcqDat.Nvd = 1;
    end
        
    S_array = reshape(S,AcqDat.Nb,AcqDat.Ndir,AcqDat.Nvd);
    nvd = 1;
    S_array = squeeze(S_array(:,:,nvd));
    Npoints = numel(S_array);
    S_vector = reshape(S_array,Npoints,1);

    indx_all_array = repmat((1:Npoints)',[1 Npoints]);
    S_array = reshape(S_vector,AcqDat.Nb,AcqDat.Ndir);
    S_PA = sum(S_array,2); S0 = max(S_PA);

    PAweight = zeros(size(S_vector));
    PAweight(td1start:end) = 1;
    PAweight_array = reshape(PAweight,AcqDat.Nb,AcqDat.Ndir);

    S_PA = sum(S_array.*PAweight_array,2)./sum(PAweight_array,2);
    S0 = max(S_PA);
    %figure(1), clf, plot(1:AcqDat.Nb,S_PA,'o'), return

    %figure(1), clf, plot(1:AcqDat.Nb,S_PA/S0,'o'), return

    Xin = AcqDat.bmat.b(1:AcqDat.Nb);
    Xin2D_iso = AcqDat.bmat.iso(1:AcqDat.Nb);
    Xin2D_aniso = AcqDat.bmat.aniso(1:AcqDat.Nb);
    
    Nb_iso = sqrt(AcqDat.Nb);
    Nb_aniso = sqrt(AcqDat.Nb);
    indx_iso = 1:Nb_iso;
    indx_aniso = 1 + Nb_iso*(0:(Nb_aniso-1));
    Xin_iso = Xin(indx_iso);
    Xin_aniso = Xin(indx_aniso);
%     Xin_iso = Xin2D_iso(indx_iso);
%     Xin_aniso = Xin2D_aniso(indx_aniso);
    %figure(1), clf, semilogy(Xin_iso,S_PA(indx_iso),'o',Xin_aniso,S_PA(indx_aniso),'o'), return

    Yin = S_PA/S0; 
    Y_iso = Yin(indx_iso);
    Y_aniso = Yin(indx_aniso);
    
    PDiso_smooth = zeros(ND_smooth,NBS);
    PDaniso_smooth = zeros(ND_smooth,NBS);
    PD_2D_smooth = zeros(NDtot_smooth,NBS);
    D_smooth = logspace(log10(Dmin),log10(Dmax),ND_smooth)';
    [Diso_smooth_a,Daniso_smooth_a] = ndgrid(D_smooth,D_smooth);
    Diso_smooth_v = reshape(Diso_smooth_a,NDtot_smooth,1);
    Daniso_smooth_v = reshape(Daniso_smooth_a,NDtot_smooth,1);

    ND_v = ceil((NDmax-NDmin)*rand(NBS,1) + NDmin);
    chisq_iso_v = zeros(NBS,1);
    chisq_aniso_v = zeros(NBS,1);
    chisq_2D_v = zeros(NBS,1);

    tic
    parfor nBS = 1:NBS
%    for nBS = 1:NBS
        nBS
    %%    
        %indx_BS = (1:Npoints)';
        indx_BS = td1start - 1 + sort(ceil((Npoints-td1start+1)*rand(Npoints,1)));
        %figure(1), clf, plot(indx_all,indx_BS,'-'), return
        indx_BS_array = repmat(indx_BS',[Npoints 1]);

        PAweight = 1e-3 + sum(indx_all_array == indx_BS_array,2);
        %figure(1), clf, plot(indx_all,PAweight,'-'), return
        PAweight_array = reshape(PAweight,AcqDat.Nb,AcqDat.Ndir);

        S_PA = sum(S_array.*PAweight_array,2)./sum(PAweight_array,2);
        %figure(1), clf, plot(1:AcqDat.Nb,S_PA,'o'), return

        Yin = S_PA/S0; 
        %figure(1), clf, plot(1:AcqDat.Nb,Yin,'o'), return
        
        Yin_iso = Yin(indx_iso);
        Yin_aniso = Yin(indx_aniso);

        ND = ND_v(nBS);
        Diso_v = 2*Dmin*(Dmax/Dmin/2^2).^rand(ND,1);
        Daniso_v = 2*Dmin*(Dmax/Dmin/2^2).^rand(ND,1);

        [K_biso1D,K_Diso1D] = ndgrid(Xin_iso,Diso_v);
        [K_baniso1D,K_Daniso1D] = ndgrid(Xin_aniso,Daniso_v);
        [K_biso2D,K_Diso2D] = ndgrid(Xin2D_iso,Diso_v);
        [K_baniso2D,K_Daniso2D] = ndgrid(Xin2D_aniso,Daniso_v);

        Kernel_iso1D = exp(-K_biso1D.*K_Diso1D);
        Kernel_aniso1D = exp(-K_baniso1D.*K_Daniso1D);
        Kernel_2D = exp(-K_biso2D.*K_Diso2D - K_baniso2D.*K_Daniso2D);

        PD_iso = lsqnonneg(Kernel_iso1D,Yin_iso);
        Ycalc_iso = Kernel_iso1D*PD_iso;
        chisq_iso = sum((Yin_iso-Ycalc_iso).^2,1)/Nb_iso;
        chisq_iso_v(nBS,1) = chisq_iso;
        PD_aniso = lsqnonneg(Kernel_aniso1D,Yin_aniso);
        Ycalc_aniso = Kernel_aniso1D*PD_aniso;
        chisq_aniso = sum((Yin_aniso-Ycalc_aniso).^2,1)/Nb_aniso;
        chisq_aniso_v(nBS,1) = chisq_aniso;
%         figure(1), clf, semilogy(Xin_iso,Yin_iso,'bo',Xin_iso,Y_iso,'bx',Xin_iso,Ycalc_iso,'b-',...
%             Xin_aniso,Yin_aniso,'ro',Xin_aniso,Y_aniso,'rx',Xin_aniso,Ycalc_aniso,'r-'), set(gca,'YLim',[1e-2 1.1])
%         pause(1)
        PD_2D = lsqnonneg(Kernel_2D,Yin);
        Ycalc_2D = Kernel_2D*PD_2D;
        chisq_2D = sum((Yin-Ycalc_2D).^2,1)/AcqDat.Nb;
        chisq_2D_v(nBS,1) = chisq_2D;

        [K_D_smooth,K_D] = ndgrid(D_smooth,Diso_v);
        Kernel_PD_smooth = exp(-(log10(K_D_smooth) - log10(K_D)).^2/2/logstd^2);
        PDtot = sum(PD_iso);
        PD_smooth = Kernel_PD_smooth*PD_iso;
        PD_iso_smooth(:,nBS) = PDtot*PD_smooth/sum(PD_smooth);

        [K_D_smooth,K_D] = ndgrid(D_smooth,Daniso_v);
        Kernel_PD_smooth = exp(-(log10(K_D_smooth) - log10(K_D)).^2/2/logstd^2);
        PDtot = sum(PD_aniso);
        PD_smooth = Kernel_PD_smooth*PD_aniso;
        PD_aniso_smooth(:,nBS) = PDtot*PD_smooth/sum(PD_smooth);
%          figure(1), clf, semilogx(D_smooth,PD_iso_smooth(:,nBS),'b-',D_smooth,PD_aniso_smooth(:,nBS),'r-'), return        
%         pause(.1)

        [K_Diso_smooth,K_Diso] = ndgrid(Diso_smooth_v,Diso_v);
        [K_Daniso_smooth,K_Daniso] = ndgrid(Daniso_smooth_v,Daniso_v);

        Kernel_PD_smooth = exp(-((log10(K_Diso_smooth) - log10(K_Diso)).^2 + (log10(K_Daniso_smooth) - log10(K_Daniso)).^2)/2/logstd^2);

        PDtot = sum(PD_2D);
        PD_smooth = Kernel_PD_smooth*PD_2D;
        PD_2D_smooth(:,nBS) = PDtot*PD_smooth/sum(PD_smooth);
%         figure(1), clf, surf(reshape(PD_smooth,sqrt(NDtot_smooth),sqrt(NDtot_smooth))), axis square, view(0,90), shading flat, colormap('hot')
%         pause(.1)

    end
    
    toc
    %%
    figure(2), clf
    semilogy(ND_v,chisq_iso_v,'bo',ND_v,chisq_aniso_v,'ro',ND_v,chisq_2D_v,'go')
    PD_iso_smooth_av = zeros(ND_smooth,NBS);
    PD_aniso_smooth_av = zeros(ND_smooth,NBS);
    PD_2D_smooth_av = zeros(NDtot_smooth,NBS);
    parfor nBS = 1:NBS
        %nBS
        indx_BS = sort(ceil(NBS*rand(NBS,1)));
        PD_iso_smooth_av(:,nBS) = mean(PD_iso_smooth(:,indx_BS),2);
        PD_aniso_smooth_av(:,nBS) = mean(PD_aniso_smooth(:,indx_BS),2);
        PD_2D_smooth_av(:,nBS) = mean(PD_2D_smooth(:,indx_BS),2);
    end

    PD_iso_smooth_std = std(PD_iso_smooth_av,0,2);
    PD_iso_smooth_mean = mean(PD_iso_smooth_av,2);

    PD_aniso_smooth_std = std(PD_aniso_smooth_av,0,2);
    PD_aniso_smooth_mean = mean(PD_aniso_smooth_av,2);

    PD_2D_smooth_std = std(PD_2D_smooth_av,0,2);
    PD_2D_smooth_mean = mean(PD_2D_smooth_av,2);

    mask = PD_iso_smooth_mean < 3*PD_iso_smooth_std;
    PD_iso_smooth_mask = PD_iso_smooth_mean;
    PD_iso_smooth_mask(mask) = 0;

    mask = PD_aniso_smooth_mean < 3*PD_aniso_smooth_std;
    PD_aniso_smooth_mask = PD_aniso_smooth_mean;
    PD_aniso_smooth_mask(mask) = 0;

    mask = PD_2D_smooth_mean < 3*PD_2D_smooth_std;
    PD_2D_smooth_mask = PD_2D_smooth_mean;
    PD_2D_smooth_mask(mask) = 0;

    [K_biso1D,K_Diso1D] = ndgrid(Xin_iso,D_smooth);
    [K_baniso1D,K_Daniso1D] = ndgrid(Xin_aniso,D_smooth);
    [K_biso2D,K_Diso2D] = ndgrid(Xin2D_iso,Diso_smooth_v);
    [K_baniso2D,K_Daniso2D] = ndgrid(Xin2D_aniso,Daniso_smooth_v);

    Kernel_iso1D = exp(-K_biso1D.*K_Diso1D);
    Kernel_aniso1D = exp(-K_baniso1D.*K_Daniso1D);
    Kernel_2D = exp(-K_biso2D.*K_Diso2D - K_baniso2D.*K_Daniso2D);
    
    PD_iso1D = PD_iso_smooth_mean;
    PD_aniso1D = PD_aniso_smooth_mean;
    PD_2D = PD_2D_smooth_mean;

    Ycalc_iso1D = Kernel_iso1D*PD_iso1D;
    Ycalc_aniso1D = Kernel_aniso1D*PD_aniso1D;
    Ycalc_2D = Kernel_2D*PD_2D;

    PD_2D = reshape(PD_2D,ND_smooth,ND_smooth);
    PD_iso2D = sum(PD_2D,2);
    PD_aniso2D = sum(PD_2D,1)';
    
    Ycalc_iso2D = Kernel_iso1D*PD_iso2D;
    Ycalc_aniso2D = Kernel_aniso1D*PD_aniso2D;
    
    Yin_iso1D = Yin(indx_iso);
    Yin_aniso1D = Yin(indx_aniso);
    
    BSdat.NBS = NBS;
    BSdat.ND = ND_smooth;
    BSdat.Diso = D_smooth;
    BSdat.Daniso = D_smooth;
    BSdat.PD_iso = PD_iso_smooth_mean;
    BSdat.PD_iso_std = PD_iso_smooth_std;
    BSdat.PD_aniso = PD_aniso_smooth_mean;
    BSdat.PD_aniso_std = PD_aniso_smooth_std;
    BSdat.PD_2D = PD_2D_smooth_mean;
    BSdat.PD_2D_std = PD_2D_smooth_std;
    BSdat.Xin_iso = Xin_iso;    
    BSdat.Xin_aniso = Xin_aniso;    
    BSdat.Yin_iso1D = Yin_iso1D;
    BSdat.Ycalc_iso1D = Ycalc_iso1D;
    BSdat.Yin_aniso1D = Yin_aniso1D;
    BSdat.Ycalc_aniso1D = Ycalc_aniso1D;
    BSdat.Yin_2D = Yin;
    BSdat.Ycalc_2D = Ycalc_2D;
    

    fs = 12;
    lw = 2;
    left = .55;
    bottom = .6;
    width = .3;
    Imin = 2e-2;

    figure(1), clf
    axh1 = axes('position',[left bottom width .35]);
    ph1 = semilogx(BSdat.Diso,BSdat.PD_iso,'b-',BSdat.Daniso,BSdat.PD_aniso,'r-');
    set(gca,'XLim',[min(BSdat.Diso) max(BSdat.Diso)],'YLim',max([BSdat.PD_iso; BSdat.PD_aniso])*[-.1 1.1] )
    set([ph1 ],'LineWidth',lw)
    axis off

    left = .12;
    axh1 = axes('position',[left bottom width .35]);
    lh1 = semilogy(Xin_iso/1e9,Ycalc_iso1D,'b-',Xin_aniso/1e9,Ycalc_aniso1D,'r-');
    hold on
    mh1 = semilogy(Xin_iso/1e9,Yin_iso1D,'bo',Xin_aniso/1e9,Yin_aniso1D,'ro');
    set([mh1],'LineWidth',.1*lw,'MarkerSize',3*lw,{'MarkerFaceColor'},{[0 0 1]; [1 0 0]})
    set([lh1],'LineWidth',.5*lw)
    set(gca,'YLim',[Imin 1.1],'XLim',max(Xin_iso/1e9)*[-.1 1.1],'TickDir','out','Box','off','LineWidth',lw,'FontSize',fs)
    xlabel('b / 10^9 sm^-^2'), ylabel('S/S_0')

    left = .12;
    bottom = .12;
    Xval = 1e-9*reshape(AcqDat.bmat.iso(1:AcqDat.Nb),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    Yval = 1e-9*reshape(AcqDat.bmat.aniso(1:AcqDat.Nb),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    axh1 = axes('position',[left bottom width .35]);
    Zval = reshape(log10(Ycalc_2D),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    Zval(Zval<log10(Imin)) = NaN;
    ph2 = plot3(Xval,Yval,Zval,'k-');
    hold on
    ph3 = plot3(Xval',Yval',Zval','k-');
    Zval = reshape(log10(Yin),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    Zval(Zval<log10(Imin)) = NaN;
    mh1 = plot3(Xval,Yval,Zval,'bo');
    view(140,10)
    axis([max(AcqDat.bmat.iso)*1e-9*[-.1 1.1 -.1 1.1] log10(Imin) .1])
    set([mh1],'LineWidth',.1*lw,'MarkerSize',3*lw,'MarkerFaceColor',[0 0 1])
    set([ph2 ph3],'LineWidth',.5*lw)
    Itick = [.01 .1 1]; ItickLabel = {'0.01', '0.1', '1'};
    set(axh1,'ZTick',log10(Itick),'ZTickLabel',ItickLabel)
    set([axh1 ],'LineWidth',lw,'FontSize',fs,'TickDir','out','Box','off')
    xlabel('b_{iso} / 10^9 sm^-^2'), ylabel('b_{aniso} / 10^9 sm^-^2'), zlabel('S/S_0')

    left = .55;
    Xval = log10(BSdat.Diso);
    Yval = log10(BSdat.Daniso);
    Zval = reshape(BSdat.PD_2D,numel(BSdat.Diso),numel(BSdat.Daniso))';
    ZprojX = sum(Zval,1);
    ZprojY = sum(Zval,2);
    axes('position',[left bottom width width*11/8.5])
    imagesc(Xval,Yval,Zval)
    axis tight, shading flat, colormap('hot')
    set(gca,'LineWidth',lw,'FontSize',fs,'YDir','normal','TickDir','out')
    xlabel('log(D_{iso} / m^2s^-^1)'), ylabel('log(D_{aniso})')
    axes('position',[left bottom+width*11/8.5+.02 width .1]);
    plot(Xval,ZprojX,'k-','LineWidth',lw)
    axis tight off
    set(gca,'YLim',max(ZprojX)*[-.02 1.02])
    axes('position',[left+width+.02 bottom .1 width*11/8.5]);
    plot(ZprojY,Yval,'k-','LineWidth',lw)
    axis tight off
    set(gca,'XLim',max(ZprojY)*[-.02 1.02])
    
    eval(['save ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/isoaniso1D2D_BSdat BSdat'])
    eval(['print -depsc -loose ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/isoaniso1D2D_BSFig'])
end

cd(wd)