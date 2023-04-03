clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/Nils_C8G2/';
ExpNam = {'C8G2_RH0_3','C8G2_RH11_3','C8G2_RH23_3','C8G2_RH33_3','C8G2_RH43_3','C8G2_RH68_3','C8G2_RH75_3','C8G2_RH96_3'};
%ExpNam = {'C8G2_RH0_3','C8G2_RH43_3','C8G2_RH75_3','C8G2_RH96_3'};
%ExpNam = {'C8G2_RH23_3'};
expno = 15:10:85;
%expno = [15 25 35 55 65 75 85];
%expno = [15 35 65 85];
ppmmin = -10; ppmmax = 130;

Peak.ppmll = [11 20 55 65 97];
Peak.ppmul = [17 38 65 86 110];
Peak.ppmll = [12 21 58 66 98];
Peak.ppmul = [16 36 64 82 107];
Peak.assignment = {'tailCH3','tail','headCH2','head','headCH'};
Peak.Ncarbons = [1 6 2 9 2];

fs = 4*11/3.33; lw = .5*1*11/3.33;

cd(DataDir)
Ndir = length(ExpNam);

Npeaks = length(Peak.ppmll);

figure(1), clf
for ndir = 1:Ndir
    ExpNam{ndir}
    cd(ExpNam{ndir})


    axes('position',[0 .15 1 .85])

    A_DP = zeros(length(expno),Npeaks);
    A_CP = zeros(length(expno),Npeaks);
    A_INEPT = zeros(length(expno),Npeaks);
    te = zeros(length(expno),1);
    
    for nexp = length(expno):-1:1
        load([num2str(expno(nexp)) '/NMRacqus.mat']);
        NMRacqus.te
        DP = load([num2str(expno(nexp)) '/SpecDat.mat']);
        CP = load([num2str(expno(nexp)+1) '/SpecDat.mat']);
        INEPT = load([num2str(expno(nexp)+2) '/SpecDat.mat']);
        
        Peak.area.DP = zeros(1,Npeaks);
        Peak.area.CP = zeros(1,Npeaks);
        Peak.area.INEPT = zeros(1,Npeaks);
        
        for npeak = 1:Npeaks
            indx = min(find(DP.Spec.ppm>Peak.ppmll(npeak))):min(find(DP.Spec.ppm>Peak.ppmul(npeak)));
            Peak.area.DP(1,npeak) = sum(DP.Spec.I(indx))/Peak.Ncarbons(npeak);
            Peak.area.CP(1,npeak) = sum(CP.Spec.I(indx))/Peak.Ncarbons(npeak);
            Peak.area.INEPT(1,npeak) = sum(INEPT.Spec.I(indx))/Peak.Ncarbons(npeak);
        end

        Imax = max([DP.Spec.I; CP.Spec.I; INEPT.Spec.I]);
        Imax = 2e6;

        figure(1)
        h = plot(DP.Spec.ppm,DP.Spec.I/Imax + (nexp-1),'k-');
        set(h,'LineWidth',1*lw)
        hold on

        h = plot(INEPT.Spec.ppm,INEPT.Spec.I/Imax + (nexp-1),'k-');
        set(h,'LineWidth',.75*lw,'Color',[1 0 0])

        h = plot(CP.Spec.ppm,CP.Spec.I/Imax + (nexp-1),'k-');
        set(h,'LineWidth',.5*lw,'Color',[0 0 1])

        h = plot([Peak.ppmll; Peak.ppmul],zeros(2,Npeaks) + (nexp-1),'-o');
        set(h,'LineWidth',.5*lw,'Color',[0 1 0])
        
        A_DP(nexp,:) = Peak.area.DP;
        A_CP(nexp,:) = Peak.area.CP;
        A_INEPT(nexp,:) = Peak.area.INEPT;
        te(nexp) = NMRacqus.te;

        eval(['save ' num2str(expno(nexp)) '/PeakDat.mat Peak']);
        %delete([num2str(expno(nexp)) '/PeakDat.mat Peak']);
    end
    figure(1)
    hold off
    set(gca,'XDir','reverse','XTick',[0 50 100],'YTick',[],...
        'LineWidth',lw,'FontSize',fs,...
        'Box','off','TickDir','out','Ycolor',[1 1 1],'Zcolor',[1 1 1])
    if ndir > 1
        axis off
        %set(gca,
    end
    
    Imin = -.1;
    Imax = length(expno)+.1;
    if exist('ppmmin') == 0
        ppmmin = min(CP.Spec.ppm);
        ppmmax = max(CP.Spec.ppm);
    end
    axis([ppmmin ppmmax Imin Imax])
    
    %%
    Peakmax = max([max(max(A_DP)) max(max(A_CP)) max(max(A_INEPT))]);
    figure(2), clf
    for npeak = 1:Npeaks
        axes('position',[.1+.9*(npeak-1)/Npeaks .2 .7*1/Npeaks .7])
        plot(te,A_DP(:,npeak)/Peakmax,'k-o')
        hold on
        plot(te,A_CP(:,npeak)/Peakmax,'b-o')
        plot(te,A_INEPT(:,npeak)/Peakmax,'r-o')
        title(Peak.assignment(npeak))
        set(gca,'YLim',[-.1 1.1])
    end

    cd ..
end

%xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'])

cd(DataDir)
%eval(['print -depsc -loose PTssNMR_MasterFig'])

cd(wd)
