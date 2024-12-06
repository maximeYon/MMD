clear all

wd = cd;

%DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
DataDir = '/Users/daniel/Dropbox/NMRdata/AVII500';
%DataDir = '/Users/daniel/Documents/Spaces/Presentations';

% ExpNam = {'RandSampDR1R2'}; expno = [30 31 33 49:55 66 68:71]; expno = [141];
ExpNam = {'AOToct_temp10'}; expno = [17:10:217]; expno = 37;
%ExpNam = {'AOToct8'}; expno = 16;
%ExpNam = {'AOToct9'}; expno = 20;
%ExpNam = {'C14E5_6'};
%ExpNam = {'C14E5_7'}; %expno = 22;
%ExpNam = {'AOToct_Eq6'}; expno = 20:10:240;
%ExpNam = {'AOToct_Eq7'}; expno = 20:10:240;

fs = 16;
lw = 2;

cd(DataDir)
cd(ExpNam{1})
if exist('expno') == 0
    GetExpnos
end
for nexp = 1:length(expno)
    load([num2str(expno(nexp)) '/NMRacqus'])
    if any(strcmp(NMRacqus.pulprog,{'DT_tbppgsteT1T2'})) == 1
        load([num2str(expno(nexp)) '/PP'])
        load([num2str(expno(nexp)) '/DTR1R2spec'])

        %uDTrec.mask_chi = find(sqrt(ReconDat.chisq) < ReconDat.chilim);
        uDTrec.mask_chi = 1:ReconDat.NBS;
        uDTrec_a.w = uDTrec_a.w(:,uDTrec.mask_chi);
        uDTrec_a.iso = uDTrec_a.iso(:,uDTrec.mask_chi);
        uDTrec_a.Delta = uDTrec_a.Delta(:,uDTrec.mask_chi);
        uDTrec_a.theta = uDTrec_a.theta(:,uDTrec.mask_chi);
        uDTrec_a.phi = uDTrec_a.phi(:,uDTrec.mask_chi);
        uDTrec_a.R1 = uDTrec_a.R1(:,uDTrec.mask_chi);
        uDTrec_a.R2 = uDTrec_a.R2(:,uDTrec.mask_chi);

        uDTrec.S0 = sum(uDTrec_a.w,1);
        uDTrec.MD = sum(uDTrec_a.w.*uDTrec_a.iso,1)...
            ./uDTrec.S0;
        uDTrec.VMD = sum(uDTrec_a.w...
            .*(uDTrec_a.iso - repmat(uDTrec.MD,[ReconDat.Nnodes 1])).^2,1)...
            ./uDTrec.S0;
        uDTrec.Ki = uDTrec.VMD./uDTrec.MD.^2;
        uDTrec.VDD = sum(uDTrec_a.w...
            .*4/5.*(uDTrec_a.iso.*uDTrec_a.Delta).^2,1)...
            ./uDTrec.S0;
        uDTrec.Ka = uDTrec.VDD./uDTrec.MD.^2;
        uDTrec.MR1 = sum(uDTrec_a.w.*uDTrec_a.R1,1)...
            ./uDTrec.S0;
        uDTrec.VR1 = sum(uDTrec_a.w...
            .*(uDTrec_a.R1 - repmat(uDTrec.MR1,[ReconDat.Nnodes 1])).^2,1)...
            ./uDTrec.S0;
        uDTrec.KR1 = uDTrec.VR1./uDTrec.MR1.^2;
        uDTrec.MR2 = sum(uDTrec_a.w.*uDTrec_a.R2,1)...
            ./uDTrec.S0;
        uDTrec.VR2 = sum(uDTrec_a.w...
            .*(uDTrec_a.R2 - repmat(uDTrec.MR2,[ReconDat.Nnodes 1])).^2,1)...
            ./uDTrec.S0;
        uDTrec.KR2 = uDTrec.VR2./uDTrec.MR2.^2;

        %mask = find(uDTrec.Ki>3);
        mask = 1:ReconDat.NBS;
        uDTrec_a.w = uDTrec_a.w(:,(mask));
        uDTrec_a.iso = uDTrec_a.iso(:,(mask));
        uDTrec_a.Delta = uDTrec_a.Delta(:,(mask));
        uDTrec_a.theta = uDTrec_a.theta(:,(mask));
        uDTrec_a.phi = uDTrec_a.phi(:,(mask));
        uDTrec_a.R1 = uDTrec_a.R1(:,(mask));
        uDTrec_a.R2 = uDTrec_a.R2(:,(mask));
        uDTrec.N = numel(mask);



    %%

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


    %%

        figure(2), clf
        
        ND = 100;
        logstd = .1;
        Dmin = ReconDat.Dmin/10^(3*logstd);
        Dmax = ReconDat.Dmax*10^(3*logstd);
        R1min = ReconDat.R1min/10^(3*logstd);
        R1max = ReconDat.R1max*10^(3*logstd);
        R2min = ReconDat.R2min/10^(3*logstd);
        R2max = ReconDat.R2max*10^(3*logstd);


        Disotick = [1e-12 1e-11 1e-10 1e-9 1e-8];
        Dratiotick = [1e-2 1e-1 1 10 100];
        R1tick = [1e-1 1 10];
        R2tick = [1 1e1 1e2 1e3];

        minclevel = .1;
        maxclevel = .9;
        Nclevel = 5;
        clevels_norm = linspace(minclevel,maxclevel,Nclevel);
        clevels_norm = logspace(log10(minclevel),log10(maxclevel),Nclevel);

        Diso = logspace(log10(Dmin),log10(Dmax),ND);
        Dratio = logspace(log10(Dmin/Dmax),log10(Dmax/Dmin),ND);
        R1 = logspace(log10(R1min),log10(R1max),ND);
        R2 = logspace(log10(R2min),log10(R2max),ND);

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
        %xlabel(XLabel), ylabel(YLabel)

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
        %xlabel(XLabel), ylabel(YLabel)
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
        %xlabel(XLabel), ylabel(YLabel)

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
        %xlabel(XLabel), ylabel(YLabel)

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
        %xlabel(XLabel), ylabel(YLabel)

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
        %xlabel(XLabel), ylabel(YLabel)
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
        %xlabel(XLabel), ylabel(YLabel)

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
        %xlabel(XLabel), ylabel(YLabel)

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
        %xlabel(XLabel), ylabel(YLabel)
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
        %xlabel(XLabel), ylabel(YLabel)
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
        %xlabel(XLabel), ylabel(YLabel)
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
        %xlabel(XLabel), ylabel(YLabel)

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
        %xlabel(XLabel), ylabel(YLabel)

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
        %xlabel(XLabel), ylabel(YLabel)
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
        %xlabel(XLabel), ylabel(YLabel)
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
        %xlabel(XLabel), ylabel(YLabel)

        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',xtick,'YTick',ytick,...
            'TickLength',ticklength)

        set(gcf, 'PaperUnits','centimeters','PaperPosition', 4*[0 0 6 6],'PaperSize', 4*[6 6]); 
        eval(['print ' num2str(expno(nexp)) '/MOAC2Fig -loose -dpdf'])

%%
        figure(3), clf
        mask_ODF = find([uDTrec.ratio>10^1]);
        %mask_ODF = find([uDTrec.ratio>10^-.2 & uDTrec.ratio<10^.2]);
        %mask_ODF = find([uDTrec.ratio<10^-1]);

        ODFdiscrete.P = uDTrec.w(mask_ODF)/sum(uDTrec.w(mask_ODF));
        ODFdiscrete.x = sin(uDTrec.theta(mask_ODF)).*cos(uDTrec.phi(mask_ODF));
        ODFdiscrete.y = sin(uDTrec.theta(mask_ODF)).*sin(uDTrec.phi(mask_ODF));
        ODFdiscrete.z = cos(uDTrec.theta(mask_ODF));

        %Smooth ODF nodes
        load UDSRTriN1000
        ODFsmooth = UDSR;
        TR = TriRep(ODFsmooth.tri, ODFsmooth.x, ODFsmooth.y, ODFsmooth.z);
        Nsubdiv = 3;
        TR=SubdivideSphericalMesh(TR,Nsubdiv);
        ODFsmooth.tri = TR.Triangulation;
        ODFsmooth.x = TR.X(:,1);
        ODFsmooth.y = TR.X(:,2);
        ODFsmooth.z = TR.X(:,3);
        ODFsmooth.N = numel(ODFsmooth.x);
        ODFsmooth.theta = acos(ODFsmooth.z);
        ODFsmooth.phi = atan2(ODFsmooth.y,ODFsmooth.x);

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
        ODFsmooth.verts = repmat(ODFsmooth.P,[1 3]).*[sin(ODFsmooth.theta).*cos(ODFsmooth.phi)...
            sin(ODFsmooth.theta).*sin(ODFsmooth.phi) ...
            cos(ODFsmooth.theta)];
        ODFsmooth.c = abs([ODFsmooth.x ODFsmooth.y ODFsmooth.z]);


        bottom = 0;
        left = 0;
        width = 1;
        axh_ODFrec = axes('position',[left bottom width width*1]);
        p = patch('Faces',ODFsmooth.tri,'Vertices',ODFsmooth.verts/max(ODFsmooth.P));
        axis tight, axis square, axis equal
        view(30,30)
        xlabel('x')
%         set(p,'FaceColor','interp','FaceVertexCData',ODFsmooth.c,...
%         'EdgeColor','none','LineWidth',1)
        set(p,'FaceColor',[1 1 1],'FaceVertexCData',ODFsmooth.c,...
        'EdgeColor','k','LineWidth',1)
        %title('ODF rec')
        axis off

        set(gcf, 'PaperUnits','centimeters','PaperPosition', 4*[0 0 6 6],'PaperSize', 4*[6 6]); 
        eval(['print ' num2str(expno(nexp)) '/MOAC2FigODF -loose -dpdf'])

        %clear uDT uDTrec uDTrec_a ReconDat
        
        ODF = ODFsmooth;
        TR = triangulation(ODF.tri, ODF.verts);
        ODF.norms = vertexNormal(TR,(1:ODF.N)');
        
        cd(num2str(expno(nexp)))
        fid = fopen(['ODF_verts.txt'], 'w');
        N = size(ODF.verts, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.5f, %8.5f, %8.5f,\n';
        for n = 1:N;
            fprintf(fid,format,ODF.verts(n,:)/max(ODF.P));
        end
        fclose(fid);

        fid = fopen(['ODF_norms.txt'], 'w');
        N = size(ODF.norms, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.5f, %8.5f, %8.5f,\n';
        for n = 1:N;
            fprintf(fid,format,ODF.norms(n,:));
        end
        fclose(fid);

        fid = fopen(['ODF_c.txt'], 'w');
        N = size(ODF.c, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.5f, %8.5f, %8.5f,\n';
        for n = 1:N;
            fprintf(fid,format,ODF.c(n,:));
        end
        fclose(fid);

        fid = fopen(['ODF_tri.txt'], 'w');
        N = size(ODF.tri, 1);
        fprintf(fid, '%8.0i,\n', N);
        format = '%8.0f, %8.0f, %8.0f,\n';
        for n = 1:N;
            fprintf(fid,format,ODF.tri(n,:)-1);
        end
        fclose(fid); 
        cd ..

    end
end
