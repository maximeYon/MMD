clear allwd = cd;cd(['/opt/topspin/data/DT/nmr']);%cd(['/Users/daniel/NMRdata/AVII200/'])%ExpNam = {'C8G2dqf'}; expno = 14;ExpNam = {'Spindiff'}; expno = 39;   AutoPhase = 1;CheckBasline = 0;CheckPeaks = 0;PlotInterm = 0;thresh = .05;td1start = 1;cd(ExpNam{1})for nexp = 1:length(expno)    ConvertAcqus = 'Y';        ConvertProcs = 'N';       MakeTextfile = 'N';    if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0        ConvertAcqus = 'Y';    end    if exist([num2str(expno(nexp)) '/pdata/1/NMRprocs.mat']) == 0        ConvertProcs = 'Y';    end    res = fExpnoInfo2(ConvertAcqus,MakeTextfile,ConvertProcs,expno(nexp));       load([num2str(expno(nexp)) '/NMRacqus.mat'])    load([num2str(expno(nexp)) '/pdata/1/NMRprocs.mat'])    eval(['load ' num2str(expno(nexp)) '/NMRacqu2s'])    load([num2str(expno(nexp)) '/pdata/1/NMRproc2s.mat'])    td = 256*ceil(NMRacqus.td/256);    td1 = NMRacqu2s.td;    if exist([num2str(expno(nexp)) '/NMRacqu3s.mat']) == 2        eval(['load ' num2str(expno(nexp)) '/NMRacqu3s'])        td1 = td1*NMRacqu3s.td;    end    tdim.tot = td/2*td1;    de = NMRacqus.de;     swh = NMRacqus.sw_h;    sw = NMRacqus.sw_h/NMRacqus.bf1;    grpdly = NMRacqus.grpdly;    bf1 = NMRacqus.bf1;    dw = 1/NMRacqus.sw_h/2;    t = de*1e-6 + 2*dw*(0:(td/2-1))';    sw1 = NMRacqu2s.sw;    swh1 = NMRacqu2s.sw*NMRacqus.bf2;    dw1 = 1/swh1/2;    t1 = 2*dw1*(0:(td1/2-1));    if exist('si1') == 0, si1 = td1; end    nu1 = swh1*linspace(0, 1, si1)'; nu1 = nu1 - nu1(si1/2+1);    if exist('si') == 0, si = td/2; end    nu = swh*linspace(0, 1, si)'; nu = nu - nu(si/2+1);    ppm = sw*linspace(-1,0,si)' + NMRprocs.offset;    if exist('ppmshift') == 0, ppmshift = 0; end    ppm = ppm + ppmshift;    ppm1 = sw1*linspace(-1,0,si1)' + NMRproc2s.offset;    if exist('ppmshift1') == 0, ppmshift1 = 0; end    ppm1 = ppm1 + ppmshift1;    fid = fopen([num2str(expno(nexp)) '/ser'],'r','ieee-le');    Std1 = fread(fid,[2,tdim.tot],'long')';    Std1 = Std1(:,1) + i*Std1(:,2);    fclose(fid);    Std1 = reshape(Std1,td/2,td1);    %figure(1), clf, plot(real(Std1)), return    S = Std1(:,1);    S((td/2-10):td/2) = 0;    %figure(1), clf, plot(real(S)), return    if exist('lb') == 0, lb = 0; end    lbfun = exp(-abs(lb*pi*(t-t(round(grpdly)))));    S = lbfun.*S;    %figure(1), clf, plot(t,real(S),t,lbfun*max(real(S))), return    grpdlycount = (1:td/2)'*2/td - floor(td/4);    zeroshiftfun = exp(-i*(NMRacqus.grpdly*2*pi*grpdlycount));    zeroshiftfun = flipdim(zeroshiftfun,1);    zeroshiftfun = fftshift(zeroshiftfun,1);    S = fft(S,td/2,1);    S = zeroshiftfun.*S;    S = ifft(S,td/2,1);    %figure(1), clf, plot(real(S)), return    if si>(td/2), S = [S(1:(td/2-round(grpdly))); zeros(si-td/2,1); S((td/2-round(grpdly)+1):(td/2))]; end    %figure(1), clf, plot(real(S)), return    S(1) = 0.5*S(1);    I = fftshift(fft(S,si));    %figure(1), clf, plot(real(I)), return    if exist('basl') == 0, basl = [.02 .1 .9 .98]; end    baslpoints = [];    for nbasl = 1:2:length(basl)-1        baslpoints = [baslpoints round(basl(nbasl)*si):round(basl(nbasl+1)*si)];    end    baslpoints = baslpoints';    pivot = max(find(abs(I) == max(abs(I))))/si;    pivotpoint = round(pivot*si);    if exist('phc0') == 0, phc0 = 20; phc1 = 10; end    if AutoPhase        Pout = fminsearch('ACME_AutoPhase',[phc0 phc1],[],I,baslpoints,pivotpoint);        phc0=Pout(1);        phc1=Pout(2);    end    phcorrfun = exp(i*(phc0*pi/180.*ones(1,si)'+phc1*pi/180*((1:si)'-pivotpoint)./si));    I = I.*phcorrfun;    I = I - mean(I(baslpoints));    if CheckBasline        figure(1), clf        subplot(2,1,1)        plot((1:si)/si,real(I),pivotpoint/si,0,'ro')        set(gca,'XDir','reverse')        subplot(2,1,2)        plot((1:si)/si,real(I),pivotpoint/si,0,'ro',baslpoints/si,zeros(size(baslpoints)),'r.')        axis([0 1 min(real(I)) max(real(I))*.05])        set(gca,'XDir','reverse')        cd .., return    end    PolyCoeff = polyfit(baslpoints,I(baslpoints),1);    I = I - polyval(PolyCoeff,(1:si)');    %figure(1), clf, plot(1:si,real(I),baslpoints,zeros(size(baslpoints)),'r.'), return    [lbfun2,dummy] = ndgrid(lbfun,1:td1);    [zeroshiftfun2,dummy] = ndgrid(zeroshiftfun,1:td1);    [phcorrfun2,dummy] = ndgrid(phcorrfun,1:td1);    Std1((td/2-20):td/2,:) = 0;    Std1 = lbfun2.*Std1;    Std1 = fft(Std1,td/2,1);    Std1 = zeroshiftfun2.*Std1;    Std1 = ifft(Std1,td/2,1);    %%    ntd1 = [8:8:256];    figure(1), clf, semilogy(real(Std1(:,ntd1))), return%%    if si>(td/2), Std1 = [Std1(1:(td/2-round(grpdly)),:); zeros(si-td/2,td1); Std1((td/2-round(grpdly)+1):(td/2),:)]; end    %figure(1), clf, plot(real(S)), return    Std1(1,:) = 0.5*Std1(1,:);    Itd1 = fftshift(fft(Std1,si,1),1);    Itd1 = phcorrfun2.*Itd1;    Ibasl = mean(Itd1(baslpoints,:),1);    [dummy,Ibasl2] = ndgrid(1:si,Ibasl);    Itd1 = Itd1 - Ibasl2;    figure(1), clf, surf(real(Itd1)), shading('flat'), view(0,90), returnPeakPick           if any(strcmp(NMRacqus.pulprog,{'DT_T2se','DT_inxt1','DT_cpxt1'}) == 1)        Xdat = NMRacqus.vd;    elseif any(strcmp(NMRacqus.pulprog,{'DT_inxt1rho','DT_cpxt1rho'}) == 1)        Xdat = 1e-6*NMRacqus.vp;    end        ExpFit            fs = 20;    lw = 1;        figure(1), clf    axes('position',[0 .18 1 .5],'FontSize',fs*.8)    plot(PP.ppm,PP.Ispec,'k','LineWidth',lw)    set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,'Box','off','TickDir','out','Ycolor',[1 1 1])    axis('tight')    ylim = get(gca,'YLim');    ylim = .05*diff(ylim)*[-1 1] + ylim;    set(gca,'YLim',ylim)    xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'],'FontSize',fs)    hold on    plot(PP.peakppm(:,td1start),PP.Ipeak(:,td1start),'bo','LineWidth',lw)        width = 1/Npeaks;    height = .2;    left = 1-width;    bottom = .7;        for npeak = 1:Npeaks        axes('position',[left bottom width height],'FontSize',fs*.8)        plot(FitDat.Xin(:,npeak),FitDat.Yin(:,npeak),'bo')        hold on        plot(FitDat.Xin(:,npeak),FitDat.Yout(:,npeak),'k-','LineWidth',lw)        axis('tight')        ylim = get(gca,'YLim');        ylim = .05*diff(ylim)*[-1 1] + ylim; ylim(1) = 0;        xlim = get(gca,'XLim');        xlim = .1*diff(xlim)*[-1 1] + xlim;        set(gca,'XLim',xlim,'YLim',ylim,'YTick',[],'XTick',[],'Box','off')        title({[num2str(PP.peakppm(npeak,td1start),3) ' ppm'];...            [num2str(FitDat.R(npeak),3) ' s^-^1']},'FontSize',fs*.5)        left = left-width;    end        eval(['save ' num2str(expno(nexp)) '/RelaxDat PP FitDat'])    delete([num2str(expno(nexp)) '/ReportFig.*'])    eval(['print -depsc -loose ' num2str(expno(nexp)) '/ReportFig'])end