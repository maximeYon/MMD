clear all

wd = cd;

%cd(['/Users/daniel/NMRdata/AVII200/'])
cd(['/opt/topspin/data/DT/nmr'])
%ExpNam = {'Yeast20120919'}; expno = 37;
%ExpNam = {'Spindiff'}; %expno = 90;
ExpNam = {'Starch_DRcorr'}; expno = 83;

td1start = 1;
%tmixpoints = 1:8;

signal = 'intensity'; signal = 'area';
Imin = 1e-3;

cd(ExpNam{1})
if exist('expno') == 0
    GetExpnos
end

for nexp = 1:length(expno)
    eval(['load ' num2str(expno(nexp)) '/NMRacqus'])
    
    if any(strcmp(NMRacqus.pulprog,{'DT_dftepgse','DT_pgsefexsy','DT_dfpgsefexsy',...
            'DT_dftmpgse','DT_mdfpgse','DT_tepgse'})) == 1
        eval(['load ' num2str(expno(nexp)) '/NMRacqu2s'])
        eval(['load ' num2str(expno(nexp)) '/PeakDat'])

        tmix = NMRacqus.vd + 2*NMRacqus.d6;
        if NMRacqus.l11
            tmix = NMRacqus.vd + 2*NMRacqus.d6 + 2*(2*NMRacqus.d12 + NMRacqus.d13);
            if any(strcmp(NMRacqus.pulprog,{'DT_dftmpgse','DT_dfpgsefexsy'})) == 1
                tmix = NMRacqus.vd + 2*(2*NMRacqus.d12 + NMRacqus.d13);
            end                
        end
        if any(strcmp(NMRacqus.pulprog,{'DT_tepgse'})) == 1
            tmix = NMRacqus.vd + 2*(NMRacqus.d6 + NMRacqus.d4 + 2*NMRacqus.d2 + NMRacqus.d3);
        end
        if any(strcmp(NMRacqus.pulprog,{'DT_dftepgse'})) == 1
            tmix = NMRacqus.vd;
        end
        if any(strcmp(NMRacqus.pulprog,{'DT_dfpgsefexsy'})) == 1
            tmix = reshape([tmix'; tmix'],2*length(tmix),1);
        end
        if any(strcmp(NMRacqus.pulprog,{'DT_mdfpgse'})) == 1
            tmix = (NMRacqus.vc+.5)*(NMRacqus.d5 + 2*NMRacqus.d12 + NMRacqus.d13);
        end

        fitpoints = td1start:length(PP.X);

        clear FitDat
        FitDat.fitpoints = fitpoints;
        FitDat.signal = signal;
        if exist('tmixpoints') == 0
            FitDat.tmixpoints = 1:length(tmix);
        else
            FitDat.tmixpoints = tmixpoints;
        end

        for npeak = 1:PP.Npeaks
            Xin = PP.X(fitpoints);
            Yin = PP.Ipeak(:,npeak);
            if strcmp(signal,'area')==1
                Yin = PP.Apeak(:,npeak);
            end
            Yin = reshape(Yin,NMRacqus.l2,NMRacqu2s.td/NMRacqus.l2);
            Yin = Yin(FitDat.fitpoints,FitDat.tmixpoints);
            %figure(1), clf, semilogy(Xin,Yin,'o'), return

            PlotInterm = 'Y';

            Iarray = Yin/max(max(Yin));
            b = Xin;
            Ntmix = NMRacqu2s.td/NMRacqus.l2;
            tmix = tmix(FitDat.tmixpoints);
            %figure(1), clf, semilogy(b,Iarray,'o'), return

            Pin(1) = 1e-9;
            Pin(2) = 1e-10;
            Pin((1:length(tmix)) + 2) = 0.5*Iarray(1,:);
            Pin((1:length(tmix)) + 2 + length(tmix)) = 0.5*Iarray(1,:);

            Yin = Iarray;
            [Xin,dummy] = ndgrid(Xin,tmix);
            Funam = 'fBiexpCorr';

            Pout = Pin; Ynorm = 1*mean(mean(Yin)); Xnorm = mean(mean(Xin)); Pnorm = Pin; 
            Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,zeros(size(Pin)),[],[],tmix,Pnorm,Xnorm,Ynorm); %Matlab 6
            Ycalc = feval(Funam,Pout,Xin,tmix,ones(size(Pin)),1,1); error = Yin - Ycalc;
            Xcalc = logspace(-1+log10(min(min(Xin))),1+log10(max(max(Xin))),100)';
            Ycalc = feval(Funam,Pout,Xcalc,tmix,ones(size(Pin)),1,1);

            chisq = sum(sum(error.^2))
            %figure(1), clf, semilogy(Xin,Yin,'o',Xcalc,Ycalc,'k-'), cd .., cd .., return

            Icalcarray = Ycalc;
            Dfast = Pout(1);
            Dslow = Pout(2);
            I0fast = Pout((1:length(tmix)) + 2);
            I0slow = Pout((1:length(tmix)) + 2 + length(tmix));
            I0vector = I0slow + I0fast;

            FitDat.Xin(:,:,npeak) = Xin;
            FitDat.Yin(:,:,npeak) = Yin;
            FitDat.Xcalc(:,:,npeak) = Xcalc;
            FitDat.Ycalc(:,:,npeak) = Ycalc;
            FitDat.I0fast(1,:,npeak) = I0fast;
            FitDat.I0slow(1,:,npeak) = I0slow;
            FitDat.Dfast(1,npeak) = Dfast;
            FitDat.Dslow(1,npeak) = Dslow;

            %[dummy,I0array] = ndgrid(Xin(:,1),I0vector);
            %Inorm = Iarray./I0array;
            %[dummy,I0array] = ndgrid(Xcalc,I0vector);
            %Icalcnorm = Icalcarray./I0array;

            fs = 20;
            lw = 3;

            if PlotInterm
                figure(5), clf
                axes('position',[.15 .2 .33 .8],'FontSize',.8*fs)
                semilogy(Xcalc,Icalcarray/I0vector(1),'k-','LineWidth',lw)
                hold on
                semilogy(b,Iarray(:,1:2)/I0vector(1),'o','LineWidth',0.5*lw,'MarkerSize',8)
                semilogy(b,Iarray(:,3:Ntmix)/I0vector(1),'ro','LineWidth',0.5*lw,'MarkerSize',8)
                hold off
                %loglog(b,Inorm,'bo',Xcalc,Icalcnorm,'k-','LineWidth',linewidth)
                xlabel('\itb\rm / sm^-^2','FontSize',fs,'VerticalAlignment','top')
                %ylabel('\itI\rm/\itI\rm_0','FontSize',fontsize)
                ylabel('signal','FontSize',fs)
                set(gca,'LineWidth',lw,'Box','off','TickDir','out','TickLength',[.02 .02])
                %set(gca,'XTick',[1e7 1e8 1e9 1e10])
                %axis([.1*min(b) 2*max(b) max(max(Icalcarray))*[1e-4 1.3]])
                axis([max(b)*[-.1 1.1] max(max(Icalcarray))*[Imin 1.3]])
                axes('position',[.65 .2 .33 .4],'FontSize',.8*fs)
                semilogx(tmix,I0fast./I0vector,'ko','LineWidth',lw)
                hold on
                semilogx(tmix,I0slow./I0vector,'ro','LineWidth',.7*lw)
                hold off
                xlabel('\itt\rm_m / s','FontSize',fs)%, ylabel('\itX\rm_i_,_2','FontSize',30)
                ylabel('\itf\rm_e','FontSize',fs)
                axis([min(tmix) max(tmix) [-.1 1.1]])
                set(gca,'LineWidth',lw,'Box','off','TickDir','out','TickLength',[.02 .02])
                set(gca,'XTick',[1e-3 1e-2 0.1 1])
                
                axes('position',[.65 .6 .33 .4],'FontSize',.8*fs)
                semilogy(tmix,I0fast,'ko','LineWidth',lw)
                hold on
                semilogy(tmix,I0slow,'ro','LineWidth',.7*lw)
                hold off
                %xlabel('\itt\rm_m / s','FontSize',fs)%, ylabel('\itX\rm_i_,_2','FontSize',30)
                ylabel('\itI\rm_0','FontSize',fs)
                axis([min(tmix) max(tmix) [Imin*10 1.1]])
                set(gca,'LineWidth',lw,'Box','off','TickDir','out','TickLength',[.02 .02])
                set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1])
                pause(1)
            end


        end 

        delete([num2str(expno(nexp)) '/ReportFig.*'])
        eval(['print -depsc -loose ' num2str(expno(nexp)) '/ReportFig'])

        eval(['load ' num2str(expno(1)) '/NMRacqus'])
    %%
        fs = 20;
        lw = 1;

        figure(1), clf
        axes('position',[-.01 .15 1.02 .35],'FontSize',fs*.8)
        plot(PP.ppm,PP.Ispec,'k','LineWidth',lw)
        set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,'Box','off',...
            'TickDir','out','Ycolor',[1 1 1])
        axis('tight')
        ylim = get(gca,'YLim');
        ylim = .05*diff(ylim)*[-1 1] + ylim;
        set(gca,'YLim',ylim)
        xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'],'FontSize',fs)
        hold on
        plot(PP.peakppm(td1start,:),PP.Ipeak(td1start,:),'bo','LineWidth',lw)
        plot(PP.ppm,-1*ylim(1)+50*PP.Ispec,'k','LineWidth',.5*lw)
        plot(PP.peakppm(td1start,:),-1*ylim(1)+50*PP.Ipeak(td1start,:),'bo','LineWidth',lw)

        dleft = (1-.1)/PP.Npeaks;
        width = .9*dleft;
        height = .3;
        left = 1-width;
        bottom = .6;
        if PP.Npeaks == 1, width = .6*1/2; end

        for npeak = 1:PP.Npeaks
            axes('position',[left bottom width height],'FontSize',fs*.8)
            semilogy(FitDat.Xin(:,:,npeak),FitDat.Yin(:,:,npeak),'bo','LineWidth',lw)
            hold on
            semilogy(FitDat.Xcalc(:,:,npeak),FitDat.Ycalc(:,:,npeak),'k-','LineWidth',lw)

            axis('tight')
            xlabel('b / m^-^2s')
            ylim = get(gca,'YLim');
            %ylim = .05*abs(diff(ylim))*[-1 1] + abs(ylim); ylim(1) = 0;
            ylim = ylim(2)*[Imin 1.2];
            xlim = get(gca,'XLim');
            xlim = .05*diff(xlim)*[-1 1] + xlim;
            xlim = max(b)*[-.1 1.1];
            set(gca,'XLim',xlim,'YLim',ylim,...
                'TickLength',.02*[1 1],'TickDir','out','Box','off','YMinorTick','on','LineWidth',1.5*lw)
            set(gca,'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1])
            %set(gca,'XTick',[1e7 1e8 1e9 1e10 1e11])
%             title({[num2str(PP.peakppm(PP.td1start,npeak),3) ' ppm  D_e/D_i=' num2str(FitDat.Dfast(npeak)/FitDat.Dslow(npeak),2) ];
%                 ['P_i=' num2str(FitDat.Psloweq(npeak),2) ' k_i=' num2str(FitDat.ki(npeak),2) ' s^-^1']},...
%                 'FontSize',fs*.5)
            left = left-dleft;
            xlim = [min(PP.X)/2 max(PP.X)*2];
            %set(gca,'XScale','log','XLim',xlim)
            if npeak < PP.Npeaks
                set(gca,'XTick',[],'YTick',[])
            end
        end
    %%
        %eval(['save ' num2str(expno(nexp)) '/FitDat FitDat'])
    end
end