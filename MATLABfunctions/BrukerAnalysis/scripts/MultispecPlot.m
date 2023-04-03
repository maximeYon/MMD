clear all

wd = cd;

cd(['/Users/daniel/NMRdata/AVII500/DT/'])
ExpNam = {'CPMAS4mm_setup'};
expno = 32;
%ppmmin = 0; ppmmax = 85;scale = 6e-8;

ExpNam = {'C10E4relax'};
expno = 22;
scale = 1.5e-8;
expno = 84;


cd('/Users/daniel/NMRdata/AVII500/Nils_C8G2/')
ExpNam = {'C8G2_RH0_3','C8G2_RH11_3','C8G2_RH23_3','C8G2_RH33_3','C8G2_RH43_3','C8G2_RH68_3','C8G2_RH75_3','C8G2_RH96_3'};
expno = 11:10:81;



fontsize = 20;
linewidth = 2;

nexp = length(expno);

for ndir = 1:length(ExpNam)
    ExpNam{ndir}
    cd(ExpNam{ndir})


    figure(1), clf
    axes('position',[0 .15 1 .85],'LineWidth',linewidth,'FontSize',fontsize)

    load([num2str(expno(1)) '/NMRacqus.mat']);
    for nexp = length(expno):-1:1
        DP = load([num2str(expno(nexp)) '/SpecDat.mat']);
        CP = load([num2str(expno(nexp)+3) '/SpecDat.mat']);
        INEPT = load([num2str(expno(nexp)+7) '/SpecDat.mat']);

        Imax = max([DP.Spec.I; CP.Spec.I; INEPT.Spec.I])

        figure(1)
        h = plot(DP.Spec.ppm,DP.Spec.I/Imax + (nexp-1),'k-');
        set(h,'LineWidth',1.5*linewidth)
        hold on

        h = plot(INEPT.Spec.ppm,INEPT.Spec.I/Imax + (nexp-1),'k-');
        set(h,'LineWidth',.75*linewidth,'Color',[1 0 0])


        h = plot(CP.Spec.ppm,CP.Spec.I/Imax + (nexp-1),'k-');
        set(h,'LineWidth',.5*linewidth,'Color',[0 0 1])

        delete([num2str(expno(nexp)+0) '/ReportFig.*'])
        delete([num2str(expno(nexp)+1) '/ReportFig.*'])
        delete([num2str(expno(nexp)+2) '/ReportFig.*'])

    end
    figure(1)
    hold off
    set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*linewidth,...
        'Box','off','TickDir','out','Ycolor',[1 1 1],'Zcolor',[1 1 1])
    Imin = -.1;
    Imax = length(expno)+.1;
    if exist('ppmmin') == 0
        ppmmin = min(CP.Spec.ppm);
        ppmmax = max(CP.Spec.ppm);
    end
    axis([ppmmin ppmmax Imin Imax])
    xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'])
    eval(['print -depsc -loose ' num2str(expno(1)) '/ReportFig'])

    cd ..
end

cd(wd)
