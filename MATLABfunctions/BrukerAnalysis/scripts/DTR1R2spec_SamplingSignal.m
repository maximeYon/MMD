clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
%DataDir = '/Users/daniel/Dropbox';
%DataDir = '/Users/daniel/Documents/Spaces/Presentations';

% ExpNam = {'RandSampDR1R2'}; expno = [30 31 33 49:55 66 68:71]; expno = [141];
ExpNam = {'AOToct_temp10'}; expno = [17:10:217]; expno = [37];
%ExpNam = {'AOToct8'}; expno = 16;
%ExpNam = {'AOToct9'}; expno = 20;
%ExpNam = {'C14E5_6'};
%ExpNam = {'C14E5_7'}; %expno = 22;
%ExpNam = {'AOToct_Eq6'}; expno = 20:10:240;
%ExpNam = {'AOToct_Eq7'}; expno = 20:10:240;

lw = 1;
fs = 10;

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

        bT.log10TR = log10(bT.TR);
        bT.log10TE = log10(bT.TE);
        bT.log10trace = log10(bT.trace);

        left = .09;
        bottom = .12;
        width = .91;
        height = .1;
        dbottom = .15;  

        ParaNam = {'TR','log10TE','log10trace','Delta','theta','phi'};

        NPara = numel(ParaNam);

        left1 = 0.01;
        right1 = .12;
        bottom1 = .08;
        dwidth = (1-right1)/NPara;
        width = dwidth - 0.02;
        dheight = dwidth;
        height = width;

        figure(2), clf

        for nParaX = 1:NPara
            for nParaY = 1:NPara
                eval(['Xval = bT.' ParaNam{nParaX} ';'])
                eval(['Yval = bT.' ParaNam{nParaY} ';'])
                left = left1 + (nParaX-1)*dwidth;
                bottom = bottom1 + (nParaY-1)*dheight;

                axes('position',[left bottom width height])
                plot(Xval,Yval,'k.')
                axis tight
                xlim = get(gca,'XLim'); ylim = get(gca,'YLim');
                xlim = xlim + .3*abs(diff(xlim)+.1)*[-1 1];
                ylim = ylim + .3*abs(diff(ylim)+.1)*[-1 1];
                set(gca,'XLim',xlim,'YLim',ylim,'YAxisLocation','right')
                set(gca,'TickDir','out','TickLength',.03*[1 1],'LineWidth',lw,'FontSize',fs)
                if nParaX < NPara
                    set(gca,'YTickLabel',[])
                end
                if nParaY > 1
                    set(gca,'XTickLabel',[])
                end
            end
        end

        papersize = 2*8.3;
        set(gcf, 'PaperUnits','centimeters','PaperPosition', papersize*[0 0 1 1],'PaperSize', papersize*[1 1]); 
        eval(['print ' num2str(expno(nexp)) '/AcqProtocol.pdf -dpdf -loose'])


        left = 0.1;
        bottom1 = .05;
        width = .75;
        dheight = (1-.05-bottom1)/(NPara+3);
        height = dheight - 0.02;
        figure(3), clf

        Xval = 1:bT.N;
        for nPara = 1:NPara
            eval(['Yval = bT.' ParaNam{nPara} ';'])
            bottom = bottom1 + (nPara-1)*dheight;

            axes('position',[left bottom width height])
            plot(Xval,Yval,'k.')
            axis tight
            xlim = get(gca,'XLim'); ylim = get(gca,'YLim');
            xlim = xlim + .1*abs(diff(xlim)+.1)*[-1 1];
            ylim = ylim + .3*abs(diff(ylim)+.1)*[-1 1];
            set(gca,'XLim',xlim,'YLim',ylim)
            set(gca,'TickDir','out','TickLength',.005*[1 1],'LineWidth',lw,'FontSize',fs,'Box','off')
            if nPara > 1
                set(gca,'XTickLabel',[])
            end
        end


        uDTrec.mask_chi = 1:ReconDat.NBS;
        uDTrec_a.w = uDTrec_a.w(:,uDTrec.mask_chi);
        uDTrec_a.iso = uDTrec_a.iso(:,uDTrec.mask_chi);
        uDTrec_a.Delta = uDTrec_a.Delta(:,uDTrec.mask_chi);
        uDTrec_a.theta = uDTrec_a.theta(:,uDTrec.mask_chi);
        uDTrec_a.phi = uDTrec_a.phi(:,uDTrec.mask_chi);
        uDTrec_a.R1 = uDTrec_a.R1(:,uDTrec.mask_chi);
        uDTrec_a.R2 = uDTrec_a.R2(:,uDTrec.mask_chi);

        mask = 1:ReconDat.NBS;
        uDTrec_a.w = uDTrec_a.w(:,(mask));
        uDTrec_a.iso = uDTrec_a.iso(:,(mask));
        uDTrec_a.Delta = uDTrec_a.Delta(:,(mask));
        uDTrec_a.theta = uDTrec_a.theta(:,(mask));
        uDTrec_a.phi = uDTrec_a.phi(:,(mask));
        uDTrec_a.R1 = uDTrec_a.R1(:,(mask));
        uDTrec_a.R2 = uDTrec_a.R2(:,(mask));
        uDTrec.N = numel(mask);

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
        residual = (Signal-Scalc);
        chisq = sum((Signal-Scalc).^2,1)/bT.N;
        

        Yval = Signal/max(Scalc);
        bottom = bottom1 + (NPara+2)*dheight;

        axes('position',[left bottom width height])
        plot(Xval,Yval,'k.')
        axis tight
        xlim = get(gca,'XLim'); ylim = get(gca,'YLim');
        xlim = xlim + .1*abs(diff(xlim)+.1)*[-1 1];
        ylim = ylim + .3*abs(diff(ylim)+.1)*[-1 1];
        set(gca,'XLim',xlim,'YLim',ylim)
        set(gca,'TickDir','out','TickLength',.005*[1 1],'LineWidth',lw,'FontSize',fs,'Box','off')
        set(gca,'XTickLabel',[])

 
        Yval = Scalc/max(Scalc);
        bottom = bottom1 + (NPara+1)*dheight;

        axes('position',[left bottom width height])
        plot(Xval,Yval,'k.')
        axis tight
        set(gca,'XLim',xlim,'YLim',ylim)
        set(gca,'TickDir','out','TickLength',.005*[1 1],'LineWidth',lw,'FontSize',fs,'Box','off')
        set(gca,'XTickLabel',[])

        Yval = residual/max(Scalc);
        bottom = bottom1 + (NPara+0)*dheight;

        axes('position',[left bottom width height])
        plot(Xval,Yval,'k.')
        axis tight
        set(gca,'XLim',xlim,'YLim',ylim)
        set(gca,'TickDir','out','TickLength',.005*[1 1],'LineWidth',lw,'FontSize',fs,'Box','off')
        set(gca,'XTickLabel',[])

        papersize = 2*8.3;
        set(gcf, 'PaperUnits','centimeters','PaperPosition', papersize*[0 0 1 2],'PaperSize', papersize*[1 2]); 
        eval(['print ' num2str(expno(nexp)) '/AcqProtocolSignal.pdf -dpdf -loose'])

        %clear uDT uDTrec uDTrec_a ReconDat
    end
end
