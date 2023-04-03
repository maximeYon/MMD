clear all

wd = cd;

%DataDir = '/Users/daniel/Dropbox/NMRdata/AVII500';
DataDir = '/opt/topspin2/data/DT/nmr';

%ExpNam = {'Yeast_fexsy_31P'}; expno = 115;
%ExpNam = {'Yeast_sedimenation_fexsy_170302'}; expno = 2:2700; expno = 2;
%ExpNam = {'CK_170306'}; expno = 94;
%ExpNam = {'Sofia_FEXSY'}; expno = 20:419; expno = 419;
ExpNam = {'Sofia_FEXSY_20171219'}; expno = 24;

td1start = 1;
%tmixpoints = 1:8;

signal = 'area'; %signal = 'intensity'; 
Imin = 1e-3;

PlotInterm = 0;
fs = 12; lw = 1; ms = 3;

cd(DataDir)
cd(ExpNam{1})
if exist('expno') == 0
    GetExpnos
end

for nexp = 1:length(expno)
    data_path = fullfile(DataDir,ExpNam{1},num2str(expno(nexp)));
    load(fullfile(data_path,'NMRacqus'))

    if any(strcmp(NMRacqus.pulprog,{'DT_tbppgsteT1T2TM'})) == 1
        load(fullfile(data_path,'NMRacqu2s'))
        load(fullfile(data_path,'PeakDat'))

        fid = fopen(fullfile(data_path,'vdT1')); vdT1 = fscanf(fid,'%f'); fclose(fid);
        fid = fopen(fullfile(data_path,'vdT2')); vdT2 = fscanf(fid,'%f'); fclose(fid);
        fid = fopen(fullfile(data_path,'vdTM')); vdTM = fscanf(fid,'%f'); fclose(fid);

        gnams = {'ax','ay','az','bx','by','bz','cx','cy','cz'};
        for ngnam = 1:numel(gnams)
            gnam = gnams{ngnam};
            fn = fullfile(data_path,['g' gnam]);
            g = mdm_bruker_grad_read(fn);
            G.(gnam) = g;
        end

        % Load gyromagnetic ratio gamma
        gamma = mdm_bruker_gamma(NMRacqus);

        % Load max gradient Gmax
        Gmax = mdm_bruker_maxgradient(NMRacqus);
        
        %timing variables
        epsilon = NMRacqus.d2;
        tau = 2*NMRacqus.d4 + 1e-6*NMRacqus.p2;
        delta = 2*(NMRacqus.d2 + NMRacqus.d3);
        Delta = 4*NMRacqus.d2 + 2*NMRacqus.d3 + 2*NMRacqus.d4 + 2*NMRacqus.d6 + NMRacqus.d5;            
        tT2 = 12*(NMRacqus.d6 + 2*NMRacqus.d2 + NMRacqus.d3 + NMRacqus.d4) + vdT2;
        tT1 = 2*NMRacqus.d5 + vdTM;
        if NMRacqus.l11
            Delta = Delta + NMRacqus.d46 + 2*NMRacqus.d42 + NMRacqus.d43;
            tT1 = tT1 + 3*(NMRacqus.d46 + 2*NMRacqus.d42 + NMRacqus.d43);
        end
        tdiff = Delta - delta/3 - tau/2 - epsilon/2 - epsilon^2/6/delta + epsilon^3/15/delta^2;
        tmix = vdTM + 2*NMRacqus.d6;
        if NMRacqus.l11
            tmix = tmix + NMRacqus.d46 + 2*NMRacqus.d42 + NMRacqus.d43 + 2*NMRacqus.d22;
        end
        Delta_mix = tmix + 4*NMRacqus.d2 + 2*NMRacqus.d3 + 2*NMRacqus.d4;
        tdiff_mix = Delta_mix - delta/3 - tau/2 - epsilon/2 - epsilon^2/6/delta + epsilon^3/15/delta^2;
        
        G_filt = Gmax.*NMRacqus.cnst1/100*sqrt(G.ax.^2 + G.ay.^2 + G.az.^2);
        b_filt = gamma^2*G_filt.^2*delta^2*tdiff;
        G_mix = Gmax.*NMRacqus.cnst1/100*sqrt(G.bx.^2 + G.by.^2 + G.bz.^2);
        b_mix = gamma^2*G_mix.^2*delta^2.*tdiff_mix;
        G_det = Gmax.*NMRacqus.cnst1/100*sqrt(G.cx.^2 + G.cy.^2 + G.cz.^2);
        b_det = gamma^2*G_det.^2*delta^2*tdiff;
        
        N_filt = numel(unique(b_filt));
        N_mix = numel(unique(tmix));
        N_det = numel(unique(b_det));
        
        fitpoints = td1start:N_det;

        clear FitDat
        FitDat.fitpoints = fitpoints;
        FitDat.signal = signal;
        if exist('tmixpoints') == 0
            FitDat.tmixpoints = 1:N_mix;
        else
            FitDat.tmixpoints = tmixpoints;
        end
        
        PP.X = reshape(b_det,N_det,(N_mix+1));
        PP.X = PP.X(:,1);
        b = PP.X;        
        tmix = reshape(tmix,N_det,(N_mix+1));
        tmix = tmix(1,:);

        for npeak = 1:PP.Npeaks
            Xin = PP.X(fitpoints);
            Yin = PP.Ipeak(:,npeak);
            if strcmp(signal,'area')==1
                Yin = PP.Apeak(:,npeak);
            end
%             Yin = reshape(Yin,N_filt,N_det,N_mix);
%             Yin = [squeeze(Yin(1,FitDat.fitpoints,1))' squeeze(Yin(FitDat.filtpoints,FitDat.fitpoints,FitDat.tmixpoints))];
            Yin = reshape(Yin,N_det,N_mix+1);
            %figure(1), clf, semilogy(Xin,Yin,'o'), return


            Iarray = Yin/max(max(Yin));
%             b = Xin;
%             Ntmix = NMRacqu2s.td/NMRacqus.l2-1;
%             tmix = tmix(FitDat.tmixpoints);
%             Ntmix = length(tmix)-1;
%             %figure(1), clf, loglog(b,Iarray,'o'), return

            Pin(1) = 1e-9;
            Pin(2) = 5e-11;
            Pin(3) = .2;
            Pin(4) = .7;
            Pin(5) = 2;
            Pin((1:(N_mix+1)) + 5) = 1.01*Iarray(1,:);

            Yin = Iarray;
            [Xin,dummy] = ndgrid(Xin,tmix);
            Funam = 'fFEXSYfit';

            Pout = Pin; Ynorm = 1*mean(mean(Yin)); Xnorm = mean(mean(Xin)); Pnorm = Pin; 
            Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,zeros(size(Pin)),[],[],tmix,Pnorm,Xnorm,Ynorm); %Matlab 6
            Ycalc = feval(Funam,Pout,Xin,tmix,ones(size(Pin)),1,1); error = Yin - Ycalc;
            Xcalc = logspace(-1+log10(min(min(Xin))),.1+log10(max(max(Xin))),100)';
            Ycalc = feval(Funam,Pout,Xcalc,tmix,ones(size(Pin)),1,1);

            chisq = sum(sum(error.^2))
            %figure(1), clf, semilogy(Xin,Yin,'o',Xcalc,Ycalc,'k-'), cd .., cd .., return

            Icalcarray = Ycalc;
            Dfast = Pout(1);
            Dslow = Pout(2);
            Pfast0 = Pout(3);
            Pfasteq = Pout(4);
            R = Pout(5);
            I0vector = Pout(6:length(Pout));

            ki = R*Pfasteq;
            ke = R - ki;
            Psloweq = 1 - Pfasteq;
            taui = 1/ki;
            taue = 1/ke;
            H = ki*2.48e-6/3;

            FitDat.Xin(:,:,npeak) = Xin;
            FitDat.Yin(:,:,npeak) = Yin;
            FitDat.Xcalc(:,:,npeak) = Xcalc;
            FitDat.Ycalc(:,:,npeak) = Ycalc;
            FitDat.Y0(1,:,npeak) = I0vector;
            FitDat.Dfast(1,npeak) = Dfast;
            FitDat.Dslow(1,npeak) = Dslow;
            FitDat.Pfasteq(1,npeak) = Pfasteq;
            FitDat.Pfast0(1,npeak) = Pfast0;
            FitDat.ki(1,npeak) = ki;
            FitDat.ke(1,npeak) = ke;
            FitDat.Psloweq(1,npeak) = Psloweq;

            %[dummy,I0array] = ndgrid(Xin(:,1),I0vector);
            %Inorm = Iarray./I0array;
            %[dummy,I0array] = ndgrid(Xcalc,I0vector);
            %Icalcnorm = Icalcarray./I0array;

            Pfast = Pfasteq-(Pfasteq - Pfast0)*exp(-R*tmix);
            Pfast(1) = Pfasteq;

            tmixcalc = logspace(min([-2+log10(1/R) log10(tmix(1)) -2]),1+log10(1/R),100)';
            Pfastcalc = Pfasteq-(Pfasteq - Pfast0)*exp(-R*tmixcalc);


            if PlotInterm
                figure(5), clf
                axes('position',[.15 .2 .33 .8],'FontSize',.8*fs)
                semilogy(Xcalc,Icalcarray/I0vector(1),'k-','LineWidth',lw)
                hold on
                semilogy(b,Iarray(:,1:2)/I0vector(1),'o','LineWidth',0.5*lw,'MarkerSize',8)
                semilogy(b,Iarray(:,3:NMRacqu2s.td/NMRacqus.l2)/I0vector(1),'ro','LineWidth',0.5*lw,'MarkerSize',8)
                hold off
                %loglog(b,Inorm,'bo',Xcalc,Icalcnorm,'k-','LineWidth',linewidth)
                xlabel('\itb\rm / sm^-^2','FontSize',fs,'VerticalAlignment','top')
                %ylabel('\itI\rm/\itI\rm_0','FontSize',fontsize)
                ylabel('signal','FontSize',fs)
                set(gca,'LineWidth',lw,'Box','off','TickDir','out','TickLength',[.02 .02])
                %set(gca,'XTick',[1e7 1e8 1e9 1e10])
                %axis([.1*min(b) 2*max(b) max(max(Icalcarray))*[1e-4 1.3]])
                axis([max(b)*[-.1 1.1] max(max(Icalcarray))*[Imin 1.3]])
                axes('position',[.65 .2 .33 .8],'FontSize',.8*fs)
                semilogx(tmixcalc,Pfastcalc,'k-','LineWidth',lw)
                hold on
                semilogx([min(tmixcalc) max(tmixcalc)],Pfast(1)*[1 1],'--',tmix(2),Pfast(2),'o',tmix(3:Ntmix+1),Pfast(3:Ntmix+1),'o','LineWidth',lw,'MarkerSize',lw*3)
                %semilogx([min(tmixcalc) max(tmixcalc)],(1-Pfast(1))*[1 1],'--',tmix(2),(1-Pfast(2)),'o',tmix(3:Ntmix+1),(1-Pfast(3:Ntmix+1)),'o','LineWidth',linewidth,'MarkerSize',12)
                hold off
                xlabel('\itt\rm_m / s','FontSize',fs)%, ylabel('\itX\rm_i_,_2','FontSize',30)
                ylabel('\itf\rm_e','FontSize',fs)
                axis([min(tmixcalc) max(tmixcalc) [-.1 1.1]])
                set(gca,'LineWidth',lw,'XTick',[.01 .1 1],'Box','off','TickDir','out','TickLength',[.02 .02])
                pause(1)
            end


        end 


        eval(['load ' num2str(expno(1)) '/NMRacqus'])
    %%

        figure(1), clf
        axes('position',[-.01 .13 1.02 .35])
        h1 = plot(PP.ppm,PP.Ispec,'k','LineWidth',lw);
        set(gca,'XDir','reverse','YTick',[],'LineWidth',lw,'Box','off','TickDir','out','Ycolor',[1 1 1],'FontSize',fs*.8)
        axis('tight')
        ylim = get(gca,'YLim');
        ylim = .05*diff(ylim)*[-1 1] + ylim;
        set(gca,'YLim',ylim)
        hold on
        h2 = plot(PP.ppm,-1*ylim(1)+50*PP.Ispec,'k','LineWidth',.5*lw);
        set(h2,'Color',.75*[1 1 1])
        set(gca,'Children',flipud(get(gca,'Children')))
        plot(PP.peakppm(td1start,:),PP.Ipeak(td1start,:),'bo','LineWidth',lw,'MarkerSize',ms)
        plot(PP.ppm(PP.intlim),zeros(size(PP.intlim)),'r-','LineWidth',lw)
        plot(PP.peakppm(td1start,:),-1*ylim(1)+50*PP.Ipeak(td1start,:),'bo','LineWidth',lw)
        xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'],'FontSize',fs)

        dleft = (1-.1)/PP.Npeaks;
        width = .9*dleft;
        height = .3;
        left = 1-width;
        bottom = .55;
        if PP.Npeaks == 1, width = .6*1/2; end

        for npeak = 1:PP.Npeaks
            axes('position',[left bottom width height])
            if strcmp(signal,'intensity')==1
                semilogy(FitDat.Xin(:,:,npeak)/1e9,FitDat.Yin(:,:,npeak),'bo','LineWidth',lw,'MarkerSize',ms);
            elseif strcmp(signal,'area')==1
                semilogy(FitDat.Xin(:,:,npeak)/1e9,FitDat.Yin(:,:,npeak),'ro','LineWidth',lw,'MarkerSize',ms);
            end
            hold on
            semilogy(FitDat.Xcalc(:,:,npeak)/1e9,FitDat.Ycalc(:,:,npeak),'k-','LineWidth',lw)
            semilogy(max(b_filt)*[1; 1]/1e9,[Imin 1.2],'k--','LineWidth',.5*lw)

            axis('tight')
            ylim = get(gca,'YLim');
            %ylim = .05*abs(diff(ylim))*[-1 1] + abs(ylim); ylim(1) = 0;
            ylim = ylim(2)*[Imin 1.2];
            xlim = get(gca,'XLim');
            xlim = .05*diff(xlim)*[-1 1] + xlim;
            xlim = max(b)*[-.1 1.1]/1e9;
            set(gca,'XLim',xlim,'YLim',ylim,'FontSize',.5*fs,...
                'TickLength',.02*[1 1],'TickDir','out','Box','off','YMinorTick','on','LineWidth',lw)
            set(gca,'YTick',[1e-3 1e-2 1e-1 1])
            %set(gca,'XTick',[1e7 1e8 1e9 1e10 1e11])
            title({[num2str(PP.peakppm(PP.td1start,npeak),3) ' ppm'];
                ['D_e/D_i=' num2str(FitDat.Dfast(npeak)/FitDat.Dslow(npeak),2) ];
                ['P_i=' num2str(FitDat.Psloweq(npeak),2)];
                ['k_i=' num2str(FitDat.ki(npeak),2) ' s^-^1']},...
                'FontSize',fs*.5)
            xlabel('b / 10^9 m^-^2s','FontSize',.8*fs)
            ylabel('I / I_0', 'FontSize',.8*fs)
            left = left-dleft;
            xlim = [min(PP.X)/2 max(PP.X)*2];
            %set(gca,'XScale','log','XLim',xlim)
            if npeak < PP.Npeaks
                xlabel('')
                ylabel('')
                set(gca,'XTickLabel',[],'YTickLabel',[])
                title({[num2str(PP.peakppm(PP.td1start,npeak),3)];
                    [num2str(FitDat.Dfast(npeak)/FitDat.Dslow(npeak),2)];
                    [num2str(FitDat.Psloweq(npeak),2)]
                    [num2str(FitDat.ki(npeak),2)]},...
                    'FontSize',.5*fs)
            end
        end
    
        eval(['save ' num2str(expno(nexp)) '/FitDat FitDat'])
        delete([num2str(expno(nexp)) '/ReportFig.*'])
        aspect = 1.6;
        set(gcf, 'PaperUnits','centimeters','PaperPosition', 2*8.3*[0 0 1 1/aspect],'PaperSize', 2*8.3*[1 1/aspect]);
        eval(['print ' num2str(expno(nexp)) '/ReportFig -loose -dpdf'])
    end
end