clear all

wd = cd;

DataDir = '/Users/daniel/Dropbox/NMRdata/AVII500';
%DataDir = '/Users/daniel/Dropbox';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
cd(DataDir)

ExpNam = {'MOACRARE2D_test'};
%ExpNam = {'MOACRARE2D_160213'};
%ExpNam = {'Sofia_FEXSY'};
%ExpNam = {'DTD2'};
%ExpNam = {'DTD'};
% ExpNam = {'AOToct_Eq7'};
% ExpNam = {'AOToct_Eq6'};
% ExpNam = {'AOToct_Eq5'};
% ExpNam = {'AOToct_Eq4'};
% ExpNam = {'AOToct_Eq3'};
%ExpNam = {'AOToct_temp11'};
%ExpNam = {'AOToct_Eq1'};
%ExpNam = {'AOToct_temp10'};
%ExpNam = {'Pstarch_temp1'};
%ExpNam = {'Pstarch_temp2'};
%ExpNam = {'AOToct_temp9'};
%ExpNam = {'AOToct_temp8'};
%ExpNam = {'AOToct_temp7'};
%ExpNam = {'AOToct_temp'};
%ExpNam = {'AOToct_temp3'};
%ExpNam = {'AOToct_temp'};
%ExpNam = {'AOTisooctane'};
%ExpNam = {'Lalpha_3Dimag'};
%ExpNam = {'MLV_3Dimag'};

for ndir = 1:length(ExpNam)
    ExpNam{ndir}
    cd(ExpNam{ndir})

    GetExpnos

    maxexpno = max(expno);
    maxexpno = 1+floor(maxexpno/10)*10;

    % expno = 11:10:maxexpno;
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
    end
    
    time = second+60*(minute+60*(hour+24*(day+31*(month+12*year))));
    time = time-time(1);

    fs = 20; lw = 2;
    
    figure(1), clf
    axh1 = axes('position',[.2 .2 .75 .35],'FontSize',fs);
    ph1 = plot(expno,te,'-o');
    xlabel('expno'), ylabel('te')

    axh2 = axes('position',[.2 .6 .75 .35],'FontSize',fs);
    ph2 = plot(expno,time/60/60/24,'-o');
    ylabel('time / days')

    set([ph1 ph2],'LineWidth',lw)
    set([axh1 axh2],'LineWidth',lw,'TickDir','out','Box','off')
    print -dpdf TempvsExpno
    
    cd ..
end

cd(wd)