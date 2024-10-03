clear all

wd = cd;

%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/opt/topspin2/data/DT/nmr';
DataDir = '/Users/daniel/Dropbox';

ExpNam = 'Yeast_tripgse'; expno = [52:58 67 69 73 75 77]; expno = [73 75 77]; expno = 77;

td1start = 2;
signal = 'area';
PlotInterm = 0;
NBS = 1e3;
Dmin = 1e-12; Dmax = 1e-8;
NDmin = 30; NDmax = 50;
logV = .002; %variance for lognormal PD smoothing
ND_smooth = 1e2;

for nexp = 1:length(expno)
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/NMRacqus'])
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/FitDat'])
    %figure(1), clf, semilogy(FitDat.Xin{2},FitDat.Yin{2},'o'), return

    waterpeak = PP.Npeaks;
    PP.peakppm(PP.td1start,waterpeak)
    S = PP.Apeak(:,waterpeak);
    if isfield(AcqDat,'NDeltab')
        AcqDat.Nb = AcqDat.Nb*AcqDat.NDeltab;
    end
    if isfield(AcqDat,'Nvd')==0
        AcqDat.Nvd = 1;
    end
        
    %b_array = reshape(AcqDat.bmat.b,sqrt(AcqDat.Nb),sqrt(AcqDat.Nb),AcqDat.Ndir);

    %S_array = reshape(S,sqrt(AcqDat.Nb),sqrt(AcqDat.Nb),AcqDat.Ndir,AcqDat.Nvd);
    S_array = reshape(S,AcqDat.Nb,AcqDat.Ndir,AcqDat.Nvd);
    nvd = 1;
    S_array = squeeze(S_array(:,:,nvd));
    Npoints = numel(S_array);
    S_vector = reshape(S_array,Npoints,1);

    indx_all_array = repmat((1:Npoints)',[1 Npoints]);
    S_array = reshape(S_vector,AcqDat.Nb,AcqDat.Ndir);
    S_PA = sum(S_array,2); S0 = max(S_PA);
    %figure(1), clf, plot(1:AcqDat.Nb,S_PA/S0,'o'), return

    Xin = AcqDat.bmat.b(1:AcqDat.Nb);
    Xin2 = AcqDat.bmat.Delta(1:AcqDat.Nb);
    
    Nb_iso = sqrt(AcqDat.Nb);
    Nb_aniso = sqrt(AcqDat.Nb);
    indx_iso = 1:Nb_iso;
    Xin_iso = Xin(indx_iso);
    Xin2_iso = Xin2(indx_iso);
    indx_aniso = 1 + Nb_iso*(0:(Nb_aniso-1));
    Xin_aniso = Xin(indx_aniso);
    Xin2_aniso = Xin2(indx_aniso);
    %figure(1), clf, semilogy(Xin_iso,S_PA(indx_iso),'o',Xin_aniso,S_PA(indx_aniso),'o'), return
    
    PDiso_smooth = zeros(ND_smooth,NBS);
    D_smooth = logspace(log10(Dmin),log10(Dmax),ND_smooth);

    PDaniso_smooth = zeros(ND_smooth,NBS);

    ND_v = ceil((NDmax-NDmin)*rand(NBS,1) + NDmin);
    chisq_v = zeros(NBS,1);

    %tic
    p =  TimedProgressBar( NBS, 30, ...
    'Computing. Remaining time: ', ', Completed: ', 'Concluded in ' );
    parfor nBS = 1:NBS
%    for nBS = 1:NBS
        %nBS
    %%    
        %indx_BS = (1:Npoints)';
        indx_BS = td1start - 1 + sort(ceil((Npoints-td1start+1)*rand(Npoints,1)));
        %figure(1), clf, plot(indx_all,indx_BS,'-'), return
        indx_BS_array = repmat(indx_BS',[Npoints 1]);

        PAweight = 1e-3 + sum(indx_all_array == indx_BS_array,2);
        %figure(1), clf, plot(indx_all,PAweight,'-'), return
        PAweight_array = reshape(PAweight,AcqDat.Nb,AcqDat.Ndir);

        S_PA = sum(S_array.*PAweight_array,2)./sum(PAweight_array,2);
        %figure(1), clf, plot(1:AcqDat.Nb,S_PA,'o'), return

        Yin = S_PA/S0*AcqDat.Ndir; 
        %figure(1), clf, plot(1:AcqDat.Nb,Yin,'o'), return
        
        Yin_iso = Yin(indx_iso);
        Yin_aniso = Yin(indx_aniso);

        ND = ND_v(nBS);
        D_v = Dmin*(Dmax/Dmin).^rand(ND,1);

        [K_biso,K_D] = ndgrid(Xin_iso,D_v);
        [K_baniso,K_D] = ndgrid(Xin_aniso,D_v);

        Kernel_iso = exp(-K_biso.*K_D);
        Kernel_aniso = exp(-K_baniso.*K_D);

        PD_iso = lsqnonneg(Kernel_iso,Yin_iso);
        Ycalc_iso = Kernel_iso*PD_iso;
        chisq_iso = sum((Yin_iso-Ycalc_iso).^2,1)/Nb_iso;
        chisq_iso_v(nBS,1) = chisq_iso;
        PD_aniso = lsqnonneg(Kernel_aniso,Yin_aniso);
        Ycalc_aniso = Kernel_aniso*PD_aniso;
        chisq_aniso = sum((Yin_aniso-Ycalc_aniso).^2,1)/Nb_aniso;
        chisq_aniso_v(nBS,1) = chisq_aniso;
        %figure(1), clf, semilogy(Xin_iso,Yin_iso,'bo',Xin_iso,Ycalc_iso,'b-',Xin_aniso,Yin_aniso,'ro',Xin_aniso,Ycalc_aniso,'r-'), set(gca,'YLim',[1e-2 1.1]), return

        [K_D_smooth,K_D] = ndgrid(D_smooth,D_v);

        Kernel_PD_smooth = exp(-(log10(K_D_smooth) - log10(K_D)).^2/logV);

        PDtot = sum(PD_iso);
        PD_smooth = Kernel_PD_smooth*PD_iso;
        PD_iso_smooth(:,nBS) = PDtot*PD_smooth/sum(PD_smooth);
        PDtot = sum(PD_aniso);
        PD_smooth = Kernel_PD_smooth*PD_aniso;
        PD_aniso_smooth(:,nBS) = PDtot*PD_smooth/sum(PD_smooth);
%          figure(1), clf, semilogx(D_smooth,PD_iso_smooth(:,nBS),'b-',D_smooth,PD_aniso_smooth(:,nBS),'r-'), return        
%         pause(.1)
        p.progress;                 % Counter for progress Report

    end
    %toc
    p.stop;  

    %%
    figure(2), clf
    semilogy(ND_v,chisq_iso_v,'bo',ND_v,chisq_aniso_v,'ro')
    PD_iso_smooth_av = zeros(ND_smooth,NBS);
    PD_aniso_smooth_av = zeros(ND_smooth,NBS);
    parfor nBS = 1:NBS
        nBS
        indx_BS = sort(ceil(NBS*rand(NBS,1)));
        PD_iso_smooth_av(:,nBS) = mean(PD_iso_smooth(:,indx_BS),2);
        PD_aniso_smooth_av(:,nBS) = mean(PD_aniso_smooth(:,indx_BS),2);
    end

    PD_iso_smooth_std = std(PD_iso_smooth_av,0,2);
    PD_iso_smooth_mean = mean(PD_iso_smooth_av,2);

    PD_aniso_smooth_std = std(PD_aniso_smooth_av,0,2);
    PD_aniso_smooth_mean = mean(PD_aniso_smooth_av,2);
    %%
    mask = PD_iso_smooth_mean < 3*PD_iso_smooth_std;
    PD_iso_smooth_mask = PD_iso_smooth_mean;
    PD_iso_smooth_mask(mask) = 0;

    mask = PD_aniso_smooth_mean < 3*PD_aniso_smooth_std;
    PD_aniso_smooth_mask = PD_aniso_smooth_mean;
    PD_aniso_smooth_mask(mask) = 0;

    figure(1), clf
    subplot(2,2,1)
    semilogx(D_smooth,PD_iso_smooth_mean,'b-',D_smooth,PD_aniso_smooth_mean,'r-')
    title('mean')
    xlabel('D')
    subplot(2,2,2)
    semilogx(D_smooth,PD_iso_smooth_mask,'b-',D_smooth,PD_aniso_smooth_mask,'r-')
    title('mask')
    
    BSdat.NBS = NBS;
    BSdat.ND = ND_smooth;
    BSdat.D_v = D_smooth;
    BSdat.PD_iso = PD_iso_smooth_mean;
    BSdat.PD_iso_std = PD_iso_smooth_std;
    BSdat.PD_aniso = PD_aniso_smooth_mean;
    BSdat.PD_aniso_std = PD_aniso_smooth_std;

    eval(['save ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/isoaniso1D_BSdat BSdat'])
    eval(['print -depsc -loose ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/isoaniso1D_BSFig'])
end

cd(wd)