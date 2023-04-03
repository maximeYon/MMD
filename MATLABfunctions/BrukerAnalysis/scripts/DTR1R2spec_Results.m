clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
%DataDir = '/Users/daniel/Dropbox';
%DataDir = '/Users/daniel/Documents/Spaces/Presentations';

% ExpNam = {'RandSampDR1R2'}; expno = [30 31 33 49:55 66 68:71]; expno = [141];
ExpNam = {'AOToct_temp10'}; expno = [17:10:217]; expno = [47];
%ExpNam = {'AOToct8'}; expno = 16;
%ExpNam = {'AOToct9'}; expno = 20;
%ExpNam = {'C14E5_6'};
%ExpNam = {'C14E5_7'}; %expno = 22;
%ExpNam = {'AOToct_Eq6'}; expno = 20:10:240;
%ExpNam = {'AOToct_Eq7'}; expno = 20:10:240;

Itick = [.01 .1 1]; ItickLabel = {'0.01', '0.1', '1'};
Dtick = [1e-12 1e-11 1e-10 1e-9 1e-8]; DtickLabel = {'-12', '-11', '-10', '-9', '-8'};
C_iso = [0 0 1];
C_aniso = [1 0 0];
fs = 12;
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

        figure(1), clf
        bottom = .88;
        left = .1;
        width = .1;
        height = .1;
        axh_bT = axes('position',[left bottom width height]);
    %     plot(bT.trace/1e9,bT.Delta,'k.')
    %     xlabel('b / 10^9 sm^-^2'),ylabel('b_\Delta')
    %     axis tight
    %     set(gca,'XLim',max(bT.trace/1e9)*[-.1 1.1],'YLim',[-.7 1.2])
    %     semilogx(bT.trace,bT.Delta,'k.')
    %     xlabel('b / sm^-^2'),ylabel('b_\Delta')
        plot(log10(bT.trace),bT.Delta,'k.')
        xlabel('log(b / sm^-^2)'),ylabel('b_\Delta')
        axis tight
        set(gca,'YLim',[-.7 1.2])

        left = .2;
        axh_Schmidt = axes('position',[left bottom width height]);
        [X,Y] = fSchmidt(bT.dir.x,bT.dir.y,bT.dir.z);
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
        plot(X,Y,'k.')
        hold on
        h_latitude = plot(latitude.X,latitude.Y,'b-');
        h_longitude = plot(longitude.X,longitude.Y,'b-');
        axis tight equal off

        left = .37;
        axh_R = axes('position',[left bottom width height]);
    %     plot(bT.TR,bT.TE,'k.')
    %     axis tight
    %     set(gca,'XLim',max(bT.TR)*[-.1 1.1],'YLim',max(bT.TE)*[-.1 1.1])
    %     semilogy(bT.TR,bT.TE,'k.')
    %     xlabel('TR / s'),ylabel('TE / s')
        plot(bT.TR,log10(bT.TE),'k.')
        xlabel('TR / s'),ylabel('log(TE / s)')
        axis tight
        set(gca,'XLim',max(bT.TR)*[-.1 1.1])

        left = .5;
        axh_spec = axes('position',[left bottom width height]);
        plot(PP.ppm,PP.Ispec,'k','LineWidth',lw)
        set(gca,'XDir','reverse','YTick',[],'Ycolor',[1 1 1])
        axis('tight')
        ylim = get(gca,'YLim');
        ylim = .1*diff(ylim)*[-1 1] + ylim;
        set(gca,'YLim',ylim)
        xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'],'FontSize',fs)
        hold on
        plot(PP.peakppm(PP.td1start,PP.Npeaks),PP.Ipeak(PP.td1start,PP.Npeaks),'bo','LineWidth',lw)

        set([axh_bT; axh_R; axh_spec],'Box','off','LineWidth',lw,'FontSize',fs,...
            'TickDir','out','TickLength',.03*[1 1])
    %%
        bottom = .7;
        left = .05;
        width = .1;
        height = .1;
        axh_hist1 = axes('position',[left bottom width height]);
        hist(sqrt(ReconDat.chisq))
        axis tight
        hold on
        ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
        xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
        plot([xlim'],[0; 0],'k-')
        set(gca,'XLim',xlim,'YLim',ylim)
        xlabel('rms\chi')

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


        left = .17;
        axh_hist2 = axes('position',[left bottom width height]);
        hist(uDTrec.S0(mask))
        axis tight
        hold on
        ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
        xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
        plot([xlim'],[0; 0],'k-')
        set(gca,'XLim',xlim,'YLim',ylim)
        xlabel('S_0')

        left = .05;
        bottom = .1;
        axh_hist3 = axes('position',[left bottom width height]);
        hist(uDTrec.MD(mask) / 1e-9)
        axis tight
        hold on
        ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
        xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
        plot([xlim'],[0; 0],'k-')
        set(gca,'XLim',xlim,'YLim',ylim)
        xlabel('MD / 10^-^9 m^2s^-^1')

        left = .17;
        axh_hist4 = axes('position',[left bottom width height]);
        hist(uDTrec.Ki(mask))
        axis tight
        hold on
        ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
        xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
        plot([xlim'],[0; 0],'k-')
        set(gca,'XLim',xlim,'YLim',ylim)
        xlabel('K_i')

        left = .3;
        axh_hist5 = axes('position',[left bottom width height]);
        hist(uDTrec.Ka(mask))
        axis tight
        hold on
        ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
        xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
        plot([xlim'],[0; 0],'k-')
        set(gca,'XLim',xlim,'YLim',ylim)
        xlabel('K_a')

        left = .05;
        bottom = .3;
        axh_hist6 = axes('position',[left bottom width height]);
        hist(uDTrec.MR1(mask))
        axis tight
        hold on
        ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
        xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
        plot([xlim'],[0; 0],'k-')
        set(gca,'XLim',xlim,'YLim',ylim)
        xlabel('MR1 / s^-^1')

        left = .17;
        axh_hist7 = axes('position',[left bottom width height]);
        hist(uDTrec.KR1(mask))
        axis tight
        hold on
        ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
        xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
        plot([xlim'],[0; 0],'k-')
        set(gca,'XLim',xlim,'YLim',ylim)
        xlabel('KR1')

        left = .05;
        bottom = .5;
        axh_hist8 = axes('position',[left bottom width height]);
        hist(uDTrec.MR2(mask))
        axis tight
        hold on
        ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
        xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
        plot([xlim'],[0; 0],'k-')
        set(gca,'XLim',xlim,'YLim',ylim)
        xlabel('MR2 / s^-^1')

        left = .17;
        axh_hist9 = axes('position',[left bottom width height]);
        hist(uDTrec.KR2(mask))
        axis tight
        hold on
        ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
        xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
        plot([xlim'],[0; 0],'k-')
        set(gca,'XLim',xlim,'YLim',ylim)
        xlabel('KR2')

        set([axh_hist1; axh_hist2; axh_hist3; axh_hist4; axh_hist5; axh_hist6; axh_hist7; axh_hist8; axh_hist9],'Box','off','TickDir','out',...
            'TickLength',.03*[1 1],'YTick',[],'FontSize',fs,'LineWidth',lw)

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

        K.bT.trace = repmat(bT.trace,1,uDTrec.N);
        K.bT.Delta = repmat(bT.Delta,1,uDTrec.N);
        K.bT.theta = repmat(bT.theta,1,uDTrec.N);
        K.bT.phi = repmat(bT.phi,1,uDTrec.N);
        K.bT.TR = repmat(bT.TR,1,uDTrec.N);
        K.bT.TE = repmat(bT.TE,1,uDTrec.N);
        K.uDT.iso = repmat(uDTrec.iso',bT.N,1);
        K.uDT.Delta = repmat(uDTrec.Delta',bT.N,1);
        K.uDT.theta = repmat(uDTrec.theta',bT.N,1);
        K.uDT.phi = repmat(uDTrec.phi',bT.N,1);
        K.uDT.R1 = repmat(uDTrec.R1',bT.N,1);
        K.uDT.R2 = repmat(uDTrec.R2',bT.N,1);

        %Generate signal
        K.cosbeta = cos(K.bT.theta).*cos(K.uDT.theta) + sin(K.bT.theta).*sin(K.uDT.theta).*cos(K.bT.phi - K.uDT.phi);
        K.P2cosbeta = (3*K.cosbeta.^2 - 1)/2;

        K.Deff = K.uDT.iso.*(1 + 2*K.bT.Delta.*K.uDT.Delta.*K.P2cosbeta);
        Kernel = exp(-K.bT.trace.*K.Deff)...
                .*(1 - exp(-K.bT.TR.*K.uDT.R1)).*exp(-K.bT.TE.*K.uDT.R2);

        Scalc = Kernel*uDTrec.w;
        chisq = sum((Signal-Scalc).^2,1)/bT.N;
        
        indx_prolate = find(uDTrec.ratio>10);
        indx_oblate = find(uDTrec.ratio<0.1);
        indx_iso = find([uDTrec.ratio>0.1 & uDTrec.ratio<10]);
        
        Scalc_prolate = Kernel(:,indx_prolate)*uDTrec.w(indx_prolate,1);
        Scalc_oblate = Kernel(:,indx_oblate)*uDTrec.w(indx_oblate,1);
        Scalc_iso = Kernel(:,indx_iso)*uDTrec.w(indx_iso,1);

        bottom = .4;
        left = .35;
        width = .25;
        height = .35;
        axh_S = axes('position',[left bottom width height]);
        %semilogy(1:bT.N,Signal,'k-',1:bT.N,Scalc,'b.',[1; bT.N],mean(sqrt(ReconDat.chisq))*[1; 1],'k--')
        semilogy(1:bT.N,Signal,'k-',1:bT.N,Scalc,'k.',[1; bT.N],sqrt(chisq)*[1; 1],'k--')
        hold on
        semilogy(1:bT.N,Scalc_prolate,'r.',1:bT.N,Scalc_oblate,'g.',1:bT.N,Scalc_iso,'b.')       
        axis([bT.N*[-.05 1.05] mean(uDTrec.S0)*[.001 1.1]])
        set(gca,'Box','off','TickDir','out','TickLength',.01*[1 1],...
            'FontSize',fs,'LineWidth',lw)
        xlabel('N'), ylabel('S')

    %%

        ND = 100;
        logstd = .1;
        Dmin = ReconDat.Dmin/10^(3*logstd);
        Dmax = ReconDat.Dmax*10^(3*logstd);

        Diso = logspace(log10(Dmin),log10(Dmax),ND);
        Dratio = logspace(log10(Dmin/Dmax),log10(Dmax/Dmin),ND);
        [Diso_a,Dratio_a] = ndgrid(Diso,Dratio);
        Diso_v = reshape(Diso_a,ND^2,1);
        Dratio_v = reshape(Dratio_a,ND^2,1);

        [K1,K1u] = ndgrid(Diso_v,uDTrec.iso);
        [K2,K2u] = ndgrid(Dratio_v,uDTrec.ratio);

        Kernel = exp(-((log10(K1) - log10(K1u)).^2....
        + (log10(K2) - log10(K2u)).^2)/2/logstd^2);

        PD = Kernel*uDTrec.w;
        PDrec = sum(uDTrec.w)*PD/sum(PD);


    %%    
        minclevel = .1;
        maxclevel = .9;
        Nclevel = 5;
        clevels_norm = linspace(minclevel,maxclevel,Nclevel);
        clevels_norm = logspace(log10(minclevel),log10(maxclevel),Nclevel);

        Xval = log10(Diso);
        Yval = log10(Dratio);
        Zvalrec = reshape(PDrec/max(PDrec),numel(Xval),numel(Yval))';
        ZprojXrec = sum(Zvalrec,1); ZprojXrec = ZprojXrec/max(ZprojXrec);
        ZprojYrec = sum(Zvalrec,2);ZprojYrec = ZprojYrec/max(ZprojYrec);

        Dratiotick = [1e-2 1e-1 1 10 100];

        bottom = .1;
        left = .72;
        width = .2;
        height = width*11/8.5;
        height_proj = .08;    

        axh_2D = axes('position',[left bottom width height]);
        clevels = max(max(Zvalrec))*clevels_norm;
        contour(Xval,Yval,Zvalrec,clevels,'k','LineWidth',.5*lw);
        hold on
        xlabel('log(D_{iso}/m^2s^-^1)'), ylabel('log(D_{||}/D_\perp)')
    %     axes('position',[left bottom+width*11/8.5+.01 width height_proj]);
    %     plot(Xval,ZprojXrec,'k-','LineWidth',lw)
    %     hold on
    %     axis tight off
    %     ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)
        axes('position',[left-height_proj*8.5/11-.01 bottom height_proj*8.5/11 width*11/8.5]);
        plot(ZprojYrec,Yval,'k-','LineWidth',lw)
        hold on
        set(gca,'XDir','reverse')
        axis tight off
        xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    

        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',log10(Dtick),'XTickLabel',DtickLabel,'YTick',log10(Dratiotick),...
            'TickLength',.007/width*[1 1])

        R1min = ReconDat.R1min/10^(3*logstd);
        R1max = ReconDat.R1max*10^(3*logstd);

        Diso = logspace(log10(Dmin),log10(Dmax),ND);
        R1 = logspace(log10(R1min),log10(R1max),ND);
        [Diso_a,R1_a] = ndgrid(Diso,R1);
        Diso_v = reshape(Diso_a,ND^2,1);
        R1_v = reshape(R1_a,ND^2,1);

        [K1,K1u] = ndgrid(Diso_v,uDTrec.iso);
        [K2,K2u] = ndgrid(R1_v,uDTrec.R1);

        Kernel = exp(-((log10(K1) - log10(K1u)).^2....
        + (log10(K2) - log10(K2u)).^2)/2/logstd^2);

        PDisoR1 = Kernel*uDTrec.w;
        PDisoR1rec = sum(uDTrec.w)*PDisoR1/sum(PDisoR1);

        Xval = log10(Diso);
        Yval = log10(R1);
        Zvalrec = reshape(PDisoR1rec/max(PDisoR1rec),numel(Xval),numel(Yval))';
        ZprojXrec = sum(Zvalrec,1); ZprojXrec = ZprojXrec/max(ZprojXrec);
        ZprojYrec = sum(Zvalrec,2);ZprojYrec = ZprojYrec/max(ZprojYrec);

        R1tick = [1e-1 1 10];

        bottom = bottom + height + .01;

        axh_2D = axes('position',[left bottom width height]);
        clevels = max(max(Zvalrec))*clevels_norm;
        contour(Xval,Yval,Zvalrec,clevels,'k','LineWidth',.5*lw);
        hold on
        ylabel('log(R_{1} / s^-^1)')
    %     axes('position',[left bottom+width*11/8.5+.01 width height_proj]);
    %     plot(Xval,ZprojXrec,'k-','LineWidth',lw)
    %     hold on
    %     axis tight off
    %     ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)
        axes('position',[left-height_proj*8.5/11-.01 bottom height_proj*8.5/11 width*11/8.5]);
        plot(ZprojYrec,Yval,'k-','LineWidth',lw)
        hold on
        set(gca,'XDir','reverse')
        axis tight off
        xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    

        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',log10(Dtick),'XTickLabel',[],'YTick',log10(R1tick),...
            'TickLength',.007/width*[1 1])

        R2min = ReconDat.R2min/10^(3*logstd);
        R2max = ReconDat.R2max*10^(3*logstd);

        Diso = logspace(log10(Dmin),log10(Dmax),ND);
        R2 = logspace(log10(R2min),log10(R2max),ND);
        [Diso_a,R2_a] = ndgrid(Diso,R2);
        Diso_v = reshape(Diso_a,ND^2,1);
        R2_v = reshape(R2_a,ND^2,1);

        [K1,K1u] = ndgrid(Diso_v,uDTrec.iso);
        [K2,K2u] = ndgrid(R2_v,uDTrec.R2);

        Kernel = exp(-((log10(K1) - log10(K1u)).^2....
        + (log10(K2) - log10(K2u)).^2)/2/logstd^2);

        PDisoR2 = Kernel*uDTrec.w;
        PDisoR2rec = sum(uDTrec.w)*PDisoR2/sum(PDisoR2);


        Xval = log10(Diso);
        Yval = log10(R2);
        Zvalrec = reshape(PDisoR2rec/max(PDisoR2rec),numel(Xval),numel(Yval))';
        ZprojXrec = sum(Zvalrec,1); ZprojXrec = ZprojXrec/max(ZprojXrec);
        ZprojYrec = sum(Zvalrec,2);ZprojYrec = ZprojYrec/max(ZprojYrec);

        R2tick = [1 1e1 1e2 1e3];

        bottom = bottom + height + .01;

        axh_2D = axes('position',[left bottom width height]);
        clevels = max(max(Zvalrec))*clevels_norm;
        contour(Xval,Yval,Zvalrec,clevels,'k','LineWidth',.5*lw);
        hold on
        ylabel('log(R_{2} / s^-^1)')
        axes('position',[left bottom+width*11/8.5+.01 width height_proj]);
        plot(Xval,ZprojXrec,'k-','LineWidth',lw)
        hold on
        axis tight off
        ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)
        axes('position',[left-height_proj*8.5/11-.01 bottom height_proj*8.5/11 width*11/8.5]);
        plot(ZprojYrec,Yval,'k-','LineWidth',lw)
        hold on
        set(gca,'XDir','reverse')
        axis tight off
        xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)    

        set([axh_2D],'LineWidth',lw,'FontSize',fs,...
            'YDir','normal','YAxisLocation','right','TickDir','out','XTick',log10(Dtick),'XTickLabel',[],'YTick',log10(R2tick),...
            'TickLength',.007/width*[1 1])

        mask_ODF = find([uDTrec.ratio>10^1]);
        %mask_ODF = find([uDTrec.ratio<10^-1]);

        ODFdiscrete.P = uDTrec.w(mask_ODF)/sum(uDTrec.w(mask_ODF));
        ODFdiscrete.x = sin(uDTrec.theta(mask_ODF)).*cos(uDTrec.phi(mask_ODF));
        ODFdiscrete.y = sin(uDTrec.theta(mask_ODF)).*sin(uDTrec.phi(mask_ODF));
        ODFdiscrete.z = cos(uDTrec.theta(mask_ODF));

        %Smooth ODF nodes
        load UDSRTriN1000
        ODFsmooth = UDSR;
        TR = TriRep(ODFsmooth.tri, ODFsmooth.x, ODFsmooth.y, ODFsmooth.z);
        Nsubdiv = 1;
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


        bottom = .1;
        left = .5;
        width = .2;
        axh_ODFrec = axes('position',[left bottom width width*11/8.5]);
        p = patch('Faces',ODFsmooth.tri,'Vertices',ODFsmooth.verts/max(ODFsmooth.P));
        axis tight, axis square, axis equal
        view(30,30)
        xlabel('x')
        set(p,'FaceColor','interp','FaceVertexCData',ODFsmooth.c,...
        'EdgeColor','none','LineWidth',1)
        %title('ODF rec')
        axis off

        set(gcf, 'PaperPosition', [0 0 11 8.5],'PaperSize', [11 8.5]); 
        eval(['print ' num2str(expno(nexp)) '/ReportFig -loose -dpdf'])

        %clear uDT uDTrec uDTrec_a ReconDat
    end
end
