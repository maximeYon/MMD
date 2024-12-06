clear all

wd = cd;

DataDir = '/Users/daniel/Dropbox/NMRdata/AVII500';
%DataDir = '/Users/daniel/Dropbox';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
cd(DataDir)

ExpNam = {'Sofia_FEXSY'}; expno = (20:419)';

for ndir = 1:length(ExpNam)
    ExpNam{ndir}
    cd(ExpNam{ndir})


    maxexpno = max(expno);
    maxexpno = 1+floor(maxexpno/10)*10;
    
    eval(['load ' num2str(expno(1)) '/FitDat'])
    Npeaks = numel(FitDat.ki);

    % expno = 11:10:maxexpno;
    Pi = zeros(numel(expno),Npeaks);
    ki = zeros(numel(expno),Npeaks);
    te = zeros(size(expno));
    year = zeros(size(expno));
    month = zeros(size(expno));
    day = zeros(size(expno));
    hour = zeros(size(expno));
    minute = zeros(size(expno));
    second = zeros(size(expno));

    for nexp = 1:length(expno)
        %if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0
    %     if expno(nexp)>2461
    %         res = fExpnoInfo2('Y','N','N',expno(nexp));
    %     end
        eval(['load ' num2str(expno(nexp)) '/NMRacqus'])
        te(nexp) = NMRacqus.te;
        year(nexp) = NMRacqus.year;
        month(nexp) = NMRacqus.month;
        day(nexp) = NMRacqus.day;
        hour(nexp) = NMRacqus.hour;
        minute(nexp) = NMRacqus.minute;
        second(nexp) = NMRacqus.second;
        eval(['load ' num2str(expno(nexp)) '/FitDat'])
        Pi(nexp,:) = FitDat.Psloweq;
        ki(nexp,:) = FitDat.ki;
    end
    
    time = second+60*(minute+60*(hour+24*(day+31*(month+12*year))));
    time = time-time(1);

    fs = 20; lw = 2;
    
    figure(1), clf
    axh1 = axes('position',[.2 .2 .75 .35],'FontSize',fs);
    ph1 = plot(time/60/60,Pi,'o');
    xlabel('time / hours')
    ylabel('P_i')
    legend('glycerol','water'), legend('boxoff')
    set(gca,'YLim',[-.05,.35])

    axh2 = axes('position',[.2 .6 .75 .35],'FontSize',fs);
    ph2 = plot(time/60/60,ki,'o');
    set(axh2,'XTickLabel','')
    ylabel('k_i / s^-^1')
    set(gca,'YLim',[-.1,5.5])

    set([ph1 ph2],'LineWidth',.5*lw,'MarkerSize',1.5*lw)
    set([axh1 axh2],'LineWidth',lw,'TickDir','out','Box','off')

    aspect = 1.6;
    set(gcf, 'PaperUnits','centimeters','PaperPosition', 2*8.3*[0 0 1 1/aspect],'PaperSize', 2*8.3*[1 1/aspect]);
    eval(['print FEXSYvstime -loose -dpdf'])
    
    cd ..
end

cd(wd)