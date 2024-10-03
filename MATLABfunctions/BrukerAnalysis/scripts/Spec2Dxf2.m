clear allwd = cd;DataDir = '/Users/daniel/NMRdata/AVII500/DT/';%DataDir = '/Users/daniel/NMRdata/AVII200/';%DataDir = '/opt/topspin2/data/DT/nmr';%DataDir = '/opt/topspin/data/DT/nmr';cd(DataDir)ExpNam = {'Pstarch'}; expno = 9;basl = [.05 .3 .6 .95]; lb = 10; si = 4*1024;phc0 = -38.9860; phc1 = -40.4453;FIDeff = .8; PPlim = [52 119];ll = [1415 1554 1566 1581 1589 1608 1623 1636 2043 2060]';ul = ll+2; ll = ll-2;thresh = .1; td1start = 6;AutoPhase = 0;AutoPhaseAll = 0;CheckBasline = 0;FindPeaks = 0;CheckPeaks = 0;PlotInterm = 0;signal = 'area'; signal = 'intensity';cd(ExpNam{1})if exist('expno') == 0    GetExpnosendfor nexp = 1:length(expno)    ConvertAcqus = 'N';    ConvertProcs = 'N';    MakeTextfile = 'N';    if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0        ConvertAcqus = 'Y';    end    res = fExpnoInfo2(ConvertAcqus,MakeTextfile,ConvertProcs,expno(nexp));    eval(['load ' num2str(expno(nexp)) '/NMRacqus'])    if any(strcmp(NMRacqus.pulprog,{'DT_pgste','DT_pgse','DT_dfpgste','DT_fexsydpgste',...            'DT_T1ir','DT_fexsy3','DT_qMASdt','DT_dfvdcp'})) == 1                Spec2Dtd1%%                if AutoPhaseAll            AutoPhasetd1        end        PeakPick%         Iproj = real(sum(Itd1,2));%         figure(11), clf%         for ntd1 = 1:td1%             Iplot = Itd1(:,ntd1);%             %Iplot = Iplot/sum(Iplot);%             h1 = plot(PP.ppm,Iplot);%             set(h1,'Color',[ntd1/td1 0 1-ntd1/td1])%             %hold on%             set(gca,'XDir','reverse','XLim',[50 120])%             pause(.1)%         end%         return        fs = 20;        lw = 1;        %%        figure(1), clf        axes('position',[-.01 .15 1.02 .35],'FontSize',fs*.8)        plot(PP.ppm,PP.Ispec,'k','LineWidth',lw)        set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,'Box','off','TickDir','out','Ycolor',[1 1 1])        axis('tight')        ylim = get(gca,'YLim');        ylim = .05*diff(ylim)*[-1 1] + ylim;        set(gca,'YLim',ylim)        if exist('PPlim') == 1            set(gca,'XLim',PPlim)        end        xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'],'FontSize',fs)        hold on        plot(PP.peakppm(td1start,:),PP.Ipeak(td1start,:),'bo','LineWidth',lw)        %plot(PP.ppm,-1*ylim(1)+50*PP.Ispec,'k','LineWidth',.5*lw)        plot(PP.peakppm(td1start,:),-1*ylim(1)+50*PP.Ipeak(td1start,:),'bo','LineWidth',lw)        dleft = (1-.1)/Npeaks;        width = .9*dleft;        height = .3;        left = 1-width;        bottom = .6;        if Npeaks == 1, width = .8*1/2; end        for npeak = 1:Npeaks            axes('position',[left bottom width height],'FontSize',fs*.8)            plot(1:td1,PP.Ipeak(:,npeak),'bo')            hold on            axis('tight')            ylim = get(gca,'YLim');            %ylim = .05*abs(diff(ylim))*[-1 1] + abs(ylim); ylim(1) = 0;            ylim = ylim(2)*[0 1.1];            xlim = get(gca,'XLim');            xlim = .05*diff(xlim)*[-1 1] + xlim;            set(gca,'XLim',xlim,'YLim',ylim,...                'TickLength',.02*[1 1],'TickDir','out','Box','off','YMinorTick','on','LineWidth',1.5*lw)            %set(gca,'YTick',[1e-3 1e-2 1e-1 1])            %set(gca,'XTick',[1e7 1e8 1e9 1e10 1e11])            title({[num2str(PP.peakppm(td1start,npeak),3) ' ppm'];...                },'FontSize',fs*.5)            left = left-dleft;            %set(gca,'XLim',xlim,'XScale','log')            if npeak < Npeaks                set(gca,'XTick',[],'YTick',[])            end        end               %%        figure(11), clf        for npeak = 1:Npeaks            Iplot = PP.Ipeak(:,npeak);            Iplot = Iplot/sum(Iplot);            h1 = plot(1:td1,Iplot);            set(h1,'Color',[npeak/Npeaks 0 1-npeak/Npeaks])            hold on            %set(gca,'XDir','reverse','XLim',[50 120])            pause(.1)        end        %%        %eval(['save ' num2str(expno(nexp)) '/FitDat FitDat'])        delete([num2str(expno(nexp)) '/ReportFig.*'])        eval(['print -depsc -loose ' num2str(expno(nexp)) '/ReportFig'])    endend