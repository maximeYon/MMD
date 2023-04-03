clear all

wd = cd;

cd(['/Users/daniel/NMRdata/AVII500/DT/'])
% ExpNam = {'AOToct_temp2'};
% expno = 15:30:185;
% ExpNam = {'Pstarch_temp1'};
% expno = 18:10:128;
% ppmmin = 0; ppmmax = 10;
ExpNam = {'Pstarch_temp2'};
expno = [28 68 98];
ppmmin = 0; ppmmax = 10;

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
        load([num2str(expno(nexp)) '/SpecDat.mat']);
        Imax = max(Spec.I);

        figure(1)
        h = plot(Spec.ppm,Spec.I/Imax + (nexp-1),'k-');
        set(h,'LineWidth',1*linewidth)
        hold on

        delete([num2str(expno(nexp)+0) '/ReportFig.*'])
    end
    figure(1)
    hold off
    set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*linewidth,...
        'Box','off','TickDir','out','Ycolor',[1 1 1],'Zcolor',[1 1 1])
    Imin = -.1;
    Imax = length(expno)+.1;
    if exist('ppmmin') == 0
        ppmmin = min(Spec.ppm);
        ppmmax = max(Spec.ppm);
    end
    axis([ppmmin ppmmax Imin Imax])
    xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'])
    eval(['print -depsc -loose ' num2str(expno(1)) '/ReportFig'])

    cd ..
end

cd(wd)
