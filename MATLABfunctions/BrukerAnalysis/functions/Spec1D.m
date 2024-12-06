ConvertAcqus = 'N';        ConvertProcs = 'N';       MakeTextfile = 'N';
if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0
    ConvertAcqus = 'Y';
end
if exist([num2str(expno(nexp)) '/pdata/1/NMRprocs.mat']) == 0
    ConvertProcs = 'Y';
end
res = fExpnoInfo2(ConvertAcqus,MakeTextfile,ConvertProcs,expno(nexp));   

load([num2str(expno(nexp)) '/NMRacqus.mat'])
load([num2str(expno(nexp)) '/pdata/1/NMRprocs.mat'])

if exist('PPlist') == 0
    PPlist = {'zg','DT_hpdec','DT_cp','DT_inept'};
end

if any(strcmp(NMRacqus.pulprog,PPlist))
if exist([num2str(expno(nexp)) '/fid'],'file') == 2

    td = 256*ceil(NMRacqus.td/256);
    de = NMRacqus.de; 
    swh = NMRacqus.sw_h;
    sw = NMRacqus.sw_h/NMRacqus.bf1;
    grpdly = NMRacqus.grpdly;
    bf1 = NMRacqus.bf1;
    dw = 1/NMRacqus.sw_h/2;
    t = de*1e-6 + 2*dw*(0:(td/2-1))';

    if exist('si') == 0
        Proc.si = 2*td/2;
    else
        Proc.si = si;    
    end

    nu = swh*linspace(0, 1, Proc.si)'; nu = nu - nu(Proc.si/2+1);
    ppm = sw*linspace(-1,0,Proc.si)' + NMRprocs.offset;
    if exist('ppmshift') == 0
        Proc.ppmshift = 0;
    else
        Proc.ppmshift = ppmshift;    
    end
    ppm = ppm + Proc.ppmshift;

    fid = fopen([num2str(expno(nexp)) '/fid'],'r','ieee-le');
    S = fread(fid,[2,td/2],'long')';
    fclose(fid);
    
    if exist('FIDeff') == 0
        Proc.FIDeff = .9;
    else
        Proc.FIDeff = FIDeff;  
    end

    if exist('td1start') == 0
        Proc.td1start = 1;
    else
        Proc.td1start = td1start;  
    end

    S = S(:,1) + i*S(:,2);
    S(floor(Proc.FIDeff*td/2):td/2) = 0;
    %figure(1), clf, plot(real(S)), return

    if exist('lb') == 0
        Proc.lb = 0;
    else
        Proc.lb = lb;    
    end
    
    lbfun = exp(-abs(Proc.lb*pi*(t-t(round(grpdly)))));
    S = lbfun.*S;
    %figure(1), clf, plot(t,real(S),t,lbfun*max(real(S))), return

    grpdlycount = (1:td/2)'*2/td - floor(td/4);
    zeroshiftfun = exp(-i*(NMRacqus.grpdly*2*pi*grpdlycount));
    zeroshiftfun = flipdim(zeroshiftfun,1);
    zeroshiftfun = fftshift(zeroshiftfun,1);
    S = fft(S,td/2,1);
    S = zeroshiftfun.*S;
    S = ifft(S,td/2,1);
    %figure(1), clf, plot(real(S)), return

    if Proc.si>(td/2), S = [S(1:(td/2-round(grpdly))); zeros(Proc.si-td/2,1); S((td/2-round(grpdly)+1):(td/2))]; end
    %figure(1), clf, plot(real(S)), return

    S(1) = 0.5*S(1);
    I = fftshift(fft(S,Proc.si));
    %figure(1), clf, plot(real(I)), return

    if exist('basl') == 0
        Proc.basl = [.05 .15 .85 .95];
    else
        Proc.basl = basl;    
    end
    baslpoints = [];
    for nbasl = 1:2:length(Proc.basl)-1
        baslpoints = [baslpoints round(Proc.basl(nbasl)*Proc.si):round(Proc.basl(nbasl+1)*Proc.si)];
    end
    baslpoints = baslpoints';

    if exist('skip') == 0
        Proc.skip = [];
    else
        Proc.skip = skip;    
    end
    skippoints = [];
    for nskip = 1:2:length(Proc.skip)-1
        skippoints = [skippoints round(Proc.skip(nskip)*Proc.si):round(Proc.skip(nskip+1)*Proc.si)];
        skippoints(skippoints==0) = [];
    end
    skippoints = skippoints';

    if exist('pivot') == 0
        Proc.pivot = max(find(abs(I) == max(abs(I))))/Proc.si;
    else
        Proc.pivot = pivot;    
    end

    Proc.pivotpoint = round(Proc.pivot*Proc.si);

    if exist('phc0') == 0
        Proc.phc0 = 100; Proc.phc1 = 20;
    else
        Proc.phc0 = phc0; Proc.phc1 = phc1;    
    end

    if exist('AutoPhaseZeroth') == 0
        AutoPhaseZeroth = 0;
    end
    
    if AutoPhase
        if AutoPhaseZeroth
            Pout = fminsearch('ACME_AutoPhaseZeroth',[Proc.phc0],[],I,baslpoints,Proc.pivotpoint);
            Proc.phc0=Pout(1);
            Proc.phc1=0;
        else            
            Pout = fminsearch('ACME_AutoPhase',[Proc.phc0 Proc.phc1],[],I,baslpoints,Proc.pivotpoint);
            Proc.phc0=Pout(1);
            Proc.phc1=Pout(2);
        end
    end

    phcorrfun = exp(i*(Proc.phc0*pi/180.*ones(1,Proc.si)'+Proc.phc1*pi/180*((1:Proc.si)'-Proc.pivotpoint)./Proc.si));
    I = I.*phcorrfun;
    I = I - mean(I(baslpoints));


    if exist('Npolybasl') == 0
        Proc.Npolybasl = 2;
    else
        Proc.Npolybasl = Npolybasl;    
    end
    
    PolyCoeff = polyfit(baslpoints,I(baslpoints),Proc.Npolybasl);
    I = I - polyval(PolyCoeff,(1:Proc.si)');
    %figure(1), clf, plot(1:Proc.si,real(I),baslpoints,zeros(size(baslpoints)),'r.'), return

    if CheckBasline
        figure(1), clf
        subplot(2,1,1)
        plot((1:Proc.si)/Proc.si,real(I),Proc.pivotpoint/Proc.si,0,'ro')
        set(gca,'XDir','reverse')
        subplot(2,1,2)
        plot((1:Proc.si)/Proc.si,real(I),Proc.pivotpoint/Proc.si,0,'ro',baslpoints/Proc.si,zeros(size(baslpoints)),'r.')
        axis([0 1 min(real(I)) max(real(I))*.1])
        set(gca,'XDir','reverse')
        return
    end

    I = real(I);

    fs = 20;
    lw = 2;

    figure(1), clf
    axes('position',[-.01 .2 1.02 .8],'FontSize',fs*.8)
    h1 = plot(ppm,.1+10*I/max(I),'k-','LineWidth',.125*lw);
    set(h1,'Color',.75*[1 1 1])
    hold on
    plot(ppm,I/max(I),'k-','LineWidth',.5*lw)
    set(gca,'XDir','reverse','YTick',[],'LineWidth',lw,'Box','off','TickDir','out','Ycolor',[1 1 1])
    axis([min(ppm)-.02*sw max(ppm)+.02*sw 1*[min([min(I/max(I)) -.05]) 1.05]])
    xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'],'FontSize',fs)
    pause(.1)

    if MkReportFig
        set(gcf, 'PaperPosition', [0 0 11 8.5],'PaperSize', [11 8.5]); 
        eval(['print ' num2str(expno(nexp)) '/ReportFig -loose -dpdf'])
    end
    if SvSpecDataMat
        Spec.I = I;
        Spec.ppm = ppm;
        save([num2str(expno(nexp)) '/SpecDat.mat'],'Spec')
    end
    if SvSpecDataTxt
        dat = [ppm I];
        save([num2str(expno(nexp)) '/SpecDat.txt'],'dat','-ascii','-tabs')
    end
    
end
clear Proc
end
