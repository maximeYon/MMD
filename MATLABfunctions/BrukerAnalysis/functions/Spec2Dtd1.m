if exist([num2str(expno(nexp)) '/ser'])

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
    eval(['load ' num2str(expno(nexp)) '/NMRacqu2s'])
    load([num2str(expno(nexp)) '/pdata/1/NMRproc2s.mat'])

    td = 256*ceil(NMRacqus.td/256);
    td1 = NMRacqu2s.td;

    if exist([num2str(expno(nexp)) '/NMRacqu3s.mat']) == 2
        eval(['load ' num2str(expno(nexp)) '/NMRacqu3s'])
        td1 = td1*NMRacqu3s.td;
    end

    tdim.tot = td/2*td1;

    de = NMRacqus.de; 
    swh = NMRacqus.sw_h;
    sw = NMRacqus.sw_h/NMRacqus.bf1;
    grpdly = NMRacqus.grpdly;
    bf1 = NMRacqus.bf1;
    dw = 1/NMRacqus.sw_h/2;
    t = de*1e-6 + 2*dw*(0:(td/2-1))';


    sw1 = NMRacqu2s.sw;
    swh1 = NMRacqu2s.sw*NMRacqus.bf2;
    dw1 = 1/swh1/2;
    if any(strcmp(NMRacqus.pulprog,{'DT_edumbohetcorcp'})) == 1
        dw1 = NMRacqus.p20*NMRacqus.l8*0.577*1e-6/2;
        swh1 = 1/dw1/2;
        sw1 =swh1/NMRacqus.bf2;
    end
    t1 = 2*dw1*(0:(td1/2-1));
    
    if exist('si1') == 0
        Proc.si1 = td1; 
    else
        Proc.si1 = si1;
    end

    nu1 = swh1*linspace(0, 1, Proc.si1)'; nu1 = nu1 - nu1(round(Proc.si1/2+1));

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

    ppm1 = sw1*linspace(-1,0,Proc.si1)' + NMRproc2s.offset;
    if exist('ppmshift1') == 0
        Proc.ppmshift1 = 0;
    else
        Proc.ppmshift1 = ppmshift1;    
    end
    ppm1 = ppm1 + Proc.ppmshift1;

    fid = fopen([num2str(expno(nexp)) '/ser'],'r','ieee-le');
    if strcmp(NMRacqus.datapath(1:13),'/opt/topspin4')
        Std1 = fread(fid,[2,tdim.tot],'float64')';
    else
        Std1 = fread(fid,[2,tdim.tot],'long')';
    end
    Std1 = Std1(:,1) + 1i*Std1(:,2);
    fclose(fid);
    %figure(1), clf, plot(real(Std1)), return

    Std1 = reshape(Std1,td/2,td1);
    %figure(1), clf, plot(real(Std1)), return

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
    if any(strcmp(NMRacqus.pulprog,{'DT_T1ir','DT_seT1ir'})) == 1
        Proc.td1start = NMRacqu2s.td;
    end
    if any(strcmp(NMRacqus.pulprog,{'DT_imse3dcs','DT_imqe3dcs'})) == 1
        Proc.td1start = NMRacqus.l1/2*NMRacqus.l2*NMRacqus.l3 + NMRacqus.l2/2*NMRacqus.l3 + NMRacqus.l3/2+1;
    end

    S = Std1(:,Proc.td1start);
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
        Proc.basl = [.02 .1 .9 .98];
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
        Proc.phc0 = 100; Proc.phc1 = 10;
    else
        Proc.phc0 = phc0; Proc.phc1 = phc1;    
    end

    if exist('AutoPhase') == 0
        AutoPhase = 1;
    end

    if AutoPhase
        Pout = fminsearch('ACME_AutoPhase',[Proc.phc0 Proc.phc1],[],I,baslpoints,Proc.pivotpoint,skippoints);
        Proc.phc0=Pout(1);
        Proc.phc1=Pout(2);
    end
    Proc

    phcorrfun = exp(1i*(Proc.phc0*pi/180.*ones(1,Proc.si)'+Proc.phc1*pi/180*((1:Proc.si)'-Proc.pivotpoint)./Proc.si));
    I = I.*phcorrfun;
    I = I - mean(I(baslpoints));

    if CheckBasline
        figure(1), clf
        subplot(2,1,1)
        plot((1:Proc.si)/Proc.si,real(I),Proc.pivotpoint/Proc.si,0,'ro')
        set(gca,'XDir','reverse')
        subplot(2,1,2)
        plot((1:Proc.si)/Proc.si,real(I),Proc.pivotpoint/Proc.si,0,'ro',baslpoints/Proc.si,zeros(size(baslpoints)),'r.')
        axis([0 1 min(real(I)) max(real(I))*.05])
        set(gca,'XDir','reverse')
        hold on
        plot(skippoints/Proc.si,zeros(size(skippoints)),'w.')
        cd .., return
    end

    if exist('NPolyCoeff') == 0
        Proc.NPolyCoeff = 2;
    else
        Proc.NPolyCoeff = NPolyCoeff;
    end

    
    PolyCoeff = polyfit(baslpoints,I(baslpoints),Proc.NPolyCoeff);
    I = I - polyval(PolyCoeff,(1:Proc.si)');
%     figure(1), clf, plot(1:Proc.si,real(I),baslpoints,zeros(size(baslpoints)),'r.'), return

    [lbfun2,dummy] = ndgrid(lbfun,1:td1);
    [zeroshiftfun2,dummy] = ndgrid(zeroshiftfun,1:td1);
    [phcorrfun2,dummy] = ndgrid(phcorrfun,1:td1);

    Std1(floor(Proc.FIDeff*td/2):td/2,:) = 0; %removing the occasional FID spike at the very last points
    Std1 = lbfun2.*Std1;
    Std1 = fft(Std1,td/2,1);
    Std1 = zeroshiftfun2.*Std1;
    Std1 = ifft(Std1,td/2,1);
    %figure(1), clf, plot(real(Std1)), return

    if Proc.si>(td/2), Std1 = [Std1(1:(td/2-round(grpdly)),:); zeros(Proc.si-td/2,td1); Std1((td/2-round(grpdly)+1):(td/2),:)]; end
    %figure(1), clf, plot(real(S)), return

    Std1(1,:) = 0.5*Std1(1,:);
    Itd1 = fftshift(fft(Std1,Proc.si,1),1);
    Itd1 = phcorrfun2.*Itd1;

    Ibasl = mean(Itd1(baslpoints,:),1);
    [dummy,Ibasl2] = ndgrid(1:Proc.si,Ibasl);
    Itd1 = Itd1 - Ibasl2;

    %figure(1), clf, surf(real(Itd1)), shading('flat'), view(0,90), return

    %Itd1(skippoints,:) = 0;
    Std1 = ifft(ifftshift(Itd1,1),Proc.si,1);
    Std1(1,:) = 2*Std1(1,:);

    if exist('AutoPhaseAll')==0
        AutoPhaseAll = 0;
    end
    if AutoPhaseAll
        phc0guess = Proc.phc0;
        phc1guess = Proc.phc1;
        for ntd1 = 1:td1
            Itemp(:,1) = Itd1(:,ntd1);

            Pout = fminsearch('ACME_AutoPhase',[phc0guess phc1guess],[],Itemp,baslpoints,Proc.pivotpoint,skippoints);
            Proc.phc0=Pout(1);
            Proc.phc1=Pout(2);
            phcorrfun = exp(1i*(Proc.phc0*pi/180.*ones(1,Proc.si)'+Proc.phc1*pi/180*((1:Proc.si)'-Proc.pivotpoint)./Proc.si));
            Itemp = Itemp.*phcorrfun;
            Itemp = Itemp - mean(Itemp(baslpoints));

            Itd1(:,ntd1) = Itemp(:,1);

            if PlotInterm
                figure(1), clf
                plot(1:Proc.si,real(Itemp),'k')
                ax = [0 Proc.si max(real(Itemp))*[-.1 1.1]];
                xlabel('channel')
                title(['Autophase   expno=' num2str(expno(nexp)) '  td1=' num2str(ntd1) ...
                    ' phc0=' num2str(Proc.phc0,3) ' phc1=' num2str(Proc.phc1,3)])
                set(gca,'XDir','reverse')
                pause(.1)
            end
        end
    end
    
    clear Itemp
    for ntd1 = 1:td1
        Itemp(:,1) = Itd1(:,ntd1);

        PolyCoeff = polyfit(baslpoints,Itemp(baslpoints),Proc.NPolyCoeff);
        Itemp = Itemp - polyval(PolyCoeff,(1:Proc.si)');
        
        Itd1(:,ntd1) = Itemp(:,1);

        if PlotInterm
            figure(1), clf
            plot(1:Proc.si,real(Itemp),'k')
            ax = [0 Proc.si max(real(Itemp))*[-.1 1.1]];
            xlabel('channel')
            title(['abc   expno=' num2str(expno(nexp)) '  td1=' num2str(ntd1) ...
                ' PolyCoeff=' num2str(PolyCoeff,3)])
            set(gca,'XDir','reverse')
            pause(.1)
        end
    end

    clear Proc
end