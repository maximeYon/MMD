clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/Nils_C8G2/';
ExpNam = {'C8G2_RH0_3','C8G2_RH11_3','C8G2_RH23_3','C8G2_RH33_3','C8G2_RH43_3','C8G2_RH68_3','C8G2_RH75_3','C8G2_RH96_3'};
ExpNam = {'C8G2_RH0_3','C8G2_RH43_3','C8G2_RH75_3','C8G2_RH96_3'};
%ExpNam = {'C8G2_RH43_3'};
expno = 15:10:85;
expno = [15 25 35 55 65 75 85];
expno = [15 35 65 85];
ppmmin = -10; ppmmax = 130;


fs = 4*11/3.33; lw = .5*1*11/3.33;

cd(DataDir)
Ndir = length(ExpNam);

left = .1;
dleft = (1-left)/Ndir;
width = .9*dleft;

figure(1), clf
for ndir = 1:Ndir
    ExpNam{ndir}
    cd(ExpNam{ndir})


    axes('position',[left .15 width .85])

    
    for nexp = length(expno):-1:1
        load([num2str(expno(nexp)) '/NMRacqus.mat']);
        NMRacqus.te
        DP = load([num2str(expno(nexp)) '/SpecDat.mat']);
        CP = load([num2str(expno(nexp)+1) '/SpecDat.mat']);
        INEPT = load([num2str(expno(nexp)+2) '/SpecDat.mat']);

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

        delete([num2str(expno(nexp)+0) '/ReportFig.*'])
        delete([num2str(expno(nexp)+1) '/ReportFig.*'])
        delete([num2str(expno(nexp)+2) '/ReportFig.*'])

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
    
    left = left + dleft;
    
    cd ..
end

%xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'])

cd(DataDir)
eval(['print -depsc -loose PTssNMR_MasterFig'])

cd(wd)
