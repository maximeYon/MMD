clear all

wd = cd;

%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/opt/topspin2/data/DT/nmr';
DataDir = '/Users/daniel/Dropbox';

ExpNam = 'Yeast_tripgse'; expno = [52:58 67 69 73 75 77]; expno = [73 75 77];

td1start = 2;
signal = 'area';
PlotInterm = 0;
NBS = 1e4;
Dmin = 1e-12; Dmax = 1e-8;
NDmin = 50; NDmax = 100;
Nthetamin = 30; Nthetamax = 40;
logstd = .05; %standard deviation for lognormal PD smoothing
ND_smooth = 64;
NDtot_smooth = ND_smooth^2;

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

    PDparDperp_smooth = zeros(NDtot_smooth,NBS);
    Dpar_smooth = logspace(log10(Dmin),log10(Dmax),sqrt(NDtot_smooth));
    Dperp_smooth = logspace(log10(Dmin),log10(Dmax),sqrt(NDtot_smooth));
    [Dpar_smooth_a,Dperp_smooth_a] = ndgrid(Dpar_smooth,Dperp_smooth);
    Dpar_smooth_v = reshape(Dpar_smooth_a,NDtot_smooth,1);
    Dperp_smooth_v = reshape(Dperp_smooth_a,NDtot_smooth,1);

    PDisoDratio_smooth = zeros(NDtot_smooth,NBS);
    Diso_smooth = logspace(log10(Dmin),log10(Dmax),sqrt(NDtot_smooth));
    Dratio_smooth = logspace(log10(Dmin/Dmax),log10(Dmax/Dmin),sqrt(NDtot_smooth));
    [Diso_smooth_a,Dratio_smooth_a] = ndgrid(Diso_smooth,Dratio_smooth);
    Diso_smooth_v = reshape(Diso_smooth_a,NDtot_smooth,1);
    Dratio_smooth_v = reshape(Dratio_smooth_a,NDtot_smooth,1);

    PDiDa_smooth = zeros(NDtot_smooth,NBS);
    Di_smooth = logspace(log10(Dmin),log10(Dmax),sqrt(NDtot_smooth));
    Da_smooth = logspace(log10(Dmin),log10(Dmax),sqrt(NDtot_smooth));
    [Di_smooth_a,Da_smooth_a] = ndgrid(Di_smooth,Da_smooth);
    Di_smooth_v = reshape(Di_smooth_a,NDtot_smooth,1);
    Da_smooth_v = reshape(Da_smooth_a,NDtot_smooth,1);

    ND_v = ceil((NDmax-NDmin)*rand(NBS,1) + NDmin);
    Ntheta_v = ceil((Nthetamax-Nthetamin)*rand(NBS,1) + Nthetamin);
    chisq_v = zeros(NBS,1);

    tic
    parfor nBS = 1:NBS
%    for nBS = 1:NBS
        nBS
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

        ND = ND_v(nBS);
        Dpar_v = 2*Dmin*(Dmax/Dmin/2^2).^rand(ND,1);
        Dperp_v = 2*Dmin*(Dmax/Dmin/2^2).^rand(ND,1);
        %Dpar = logspace(log10(minD),log10(maxD),sqrt(ND));
        %Dperp = logspace(log10(minD),log10(maxD),sqrt(ND));
        %[Dpar_a,Dperp_a] = ndgrid(Dpar,Dperp);
        %Dpar_v = reshape(Dpar_a,ND,1);
        %Dperp_v = reshape(Dperp_a,ND,1);
        %figure(1), clf, loglog(Dpar_v,Dperp_v,'o'), return
        Diso_v = (Dpar_v + 2*Dperp_v)/3;
        DDelta_v = (Dpar_v - Dperp_v)/3./Diso_v;
        Dratio_v = Dpar_v./Dperp_v;
        %figure(1), clf, semilogx(Diso_v,DDelta_v,'o'), return

        [K_b,K_Diso] = ndgrid(Xin,Diso_v);
        [K_bDelta,K_DDelta] = ndgrid(Xin2,DDelta_v);

        K_bDDelDel = K_b.*K_Diso.*K_bDelta.*K_DDelta;
    %max(max(K_bDDelDel))
    %min(min(K_bDDelDel))

        Kernel = exp(-K_b.*K_Diso).*exp(K_bDDelDel).*...
            sqrt(pi)/2.*real(gammainc(3*K_bDDelDel,1/2)./sqrt(3*K_bDDelDel));

        indx = K_bDDelDel == 0;
        Kernel(indx) = exp(-K_b(indx).*K_Diso(indx));
        Kernel(K_b == 0) = 1;
        Kernel(K_bDDelDel < -10) = 0;
        indx = isnan(Kernel);
        Kernel(indx) = 0;
        indx = isinf(Kernel);
        Kernel(indx) = 0;
        %figure(1), clf, plot((1:AcqDat.Nb)',Kernel,'-'), set(gca,'YLim',[-.1 1.1]), return
        %figure(1), clf, plot(K_bDDelDel,Kernel,'o'), return

        PD = lsqnonneg(Kernel,Yin);
        Ycalc = Kernel*PD;
        chisq = sum((Yin-Ycalc).^2,1)/AcqDat.Nb;
        chisq_v(nBS,1) = chisq;
        %figure(1), clf, plot((1:AcqDat.Nb)',Yin,'o',(1:AcqDat.Nb)',Ycalc,'-'), set(gca,'YLim',[-.1 1.1]), return

        [K_Dpar_smooth,K_Dpar] = ndgrid(Dpar_smooth_v,Dpar_v);
        [K_Dperp_smooth,K_Dperp] = ndgrid(Dperp_smooth_v,Dperp_v);

        Kernel_PD_smooth = exp(-((log10(K_Dpar_smooth) - log10(K_Dpar)).^2 + (log10(K_Dperp_smooth) - log10(K_Dperp)).^2)/2/logstd^2);

        PDtot = sum(PD);
        PD_smooth = Kernel_PD_smooth*PD;
        PDparDperp_smooth(:,nBS) = PDtot*PD_smooth/sum(PD_smooth);

        [K_Diso_smooth,K_Diso] = ndgrid(Diso_smooth_v,Diso_v);
        [K_Dratio_smooth,K_Dratio] = ndgrid(Dratio_smooth_v,Dratio_v);

        Kernel_PD_smooth = exp(-((log10(K_Diso_smooth) - log10(K_Diso)).^2 + (log10(K_Dratio_smooth) - log10(K_Dratio)).^2)/2/logstd^2);

        PD_smooth = Kernel_PD_smooth*PD;
        PDisoDratio_smooth(:,nBS) = PDtot*PD_smooth/sum(PD_smooth);

        Ntheta = Ntheta_v(nBS);
        theta = pi/2*linspace(0,1,Ntheta)';
        Ptheta = sin(theta); Ptheta/sum(Ptheta);

%         dtheta = abs(theta(2)-theta(1));
%         theta = theta(1:Ntheta) + dtheta/2;
% 
%         [Di_a4,Da_a4,PDparDperp_a4,theta_a4] = ndgrid(Di_smooth,Da_smooth,PD,theta);
%         [Di_a4,Da_a4,Dpar_a4,theta_a4] = ndgrid(Di_smooth,Da_smooth,Dpar_v,theta);
%         [Di_a4,Da_a4,Dperp_a4,theta_a4] = ndgrid(Di_smooth,Da_smooth,Dperp_v,theta);
% 
%         Diso_a4 = (2*Dperp_a4 + Dpar_a4)/3;
%         Datheta_a4 = Dpar_a4.*cos(theta_a4).^2 + Dperp_a4.*sin(theta_a4).^2;
%         Dpar_a4 = []; Dperp_a4 = [];
%         P_a4 = PDparDperp_a4.*sin(theta_a4).*exp(-((log10(Datheta_a4)-log10(Da_a4)).^2 + (log10(Diso_a4)-log10(Di_a4)).^2)/2/logstd^2);
%         PDparDperp_a4 = []; theta_a4 = []; Datheta_a4 = []; Da_a4 = []; Diso_a4 = []; Di_a4 = [];
% 
%         P_a2 = sum(sum(P_a4,4),3);

        [Di_a4,Dpar_a4,theta_a4] = ndgrid(Di_smooth_v,Dpar_v,theta);
        [Da_a4,Dperp_a4,Ptheta_a4] = ndgrid(Da_smooth_v,Dperp_v,Ptheta);
        Diso_a4 = (2*Dperp_a4 + Dpar_a4)/3;
        Datheta_a4 = Dpar_a4.*cos(theta_a4).^2 + Dperp_a4.*sin(theta_a4).^2;

        Kernel_a4 = Ptheta_a4.*exp(-((log10(Da_a4) - log10(Datheta_a4)).^2 + (log10(Di_a4) - log10(Diso_a4)).^2)/2/logstd^2);
        Kernel_PD_smooth = sum(Kernel_a4,3);
        
        PD_smooth = Kernel_PD_smooth*PD;
        PDiDa_smooth(:,nBS) = PDtot*PD_smooth/sum(PD_smooth);
        
%         figure(1), clf, surf(reshape(PD_smooth,sqrt(NDtot_smooth),sqrt(NDtot_smooth))), axis square, view(0,90), shading flat, colormap('hot')
%         pause(.1)


    end
    
    toc
    %%
    figure(2), clf
    semilogy(ND_v,chisq_v,'o')
    PDparDperp_smooth_av = zeros(NDtot_smooth,NBS);
    PDisoDratio_smooth_av = zeros(NDtot_smooth,NBS);
    PDiDa_smooth_av = zeros(NDtot_smooth,NBS);
    parfor nBS = 1:NBS
        %nBS
        indx_BS = sort(ceil(NBS*rand(NBS,1)));
        PDparDperp_smooth_av(:,nBS) = mean(PDparDperp_smooth(:,indx_BS),2);
        PDisoDratio_smooth_av(:,nBS) = mean(PDisoDratio_smooth(:,indx_BS),2);
        PDiDa_smooth_av(:,nBS) = mean(PDiDa_smooth(:,indx_BS),2);
    end

    PDparDperp_smooth_std = std(PDparDperp_smooth_av,0,2);
    PDparDperp_smooth_mean = mean(PDparDperp_smooth_av,2);

    PDisoDratio_smooth_std = std(PDisoDratio_smooth_av,0,2);
    PDisoDratio_smooth_mean = mean(PDisoDratio_smooth_av,2);

    PDiDa_smooth_std = std(PDiDa_smooth_av,0,2);
    PDiDa_smooth_mean = mean(PDiDa_smooth_av,2);
%%
    mask = PDparDperp_smooth_mean < 3*PDparDperp_smooth_std;
    PDDparDperp_smooth_mask = PDparDperp_smooth_mean;
    PDDparDperp_smooth_mask(mask) = 0;

    mask = PDisoDratio_smooth_mean < 3*PDisoDratio_smooth_std;
    PDisoDratio_smooth_mask = PDisoDratio_smooth_mean;
    PDisoDratio_smooth_mask(mask) = 0;

    figure(1), clf
    subplot(2,2,1)
    imagesc(log10(Dpar_smooth),log10(Dperp_smooth),reshape(PDparDperp_smooth_mean,sqrt(NDtot_smooth),sqrt(NDtot_smooth))')
    axis square tight, shading flat, colormap('hot'), set(gca,'YDir','normal')
    xlabel('Dpar'), ylabel('Dperp')
    subplot(2,2,2)
    imagesc(log10(Diso_smooth),log10(Dratio_smooth),reshape(PDisoDratio_smooth_mean,sqrt(NDtot_smooth),sqrt(NDtot_smooth))')
    axis square tight, shading flat, colormap('hot'), set(gca,'YDir','normal')
    xlabel('Diso'), ylabel('Dpar/Dperp')
    subplot(2,2,3)
    imagesc(log10(Di_smooth),log10(Da_smooth),reshape(PDiDa_smooth_mean,sqrt(NDtot_smooth),sqrt(NDtot_smooth))')
    axis square tight, shading flat, colormap('hot'), set(gca,'YDir','normal')
    xlabel('Di'), ylabel('Da')

    BSdat.NBS = NBS;
    BSdat.ND = NDtot_smooth;
    BSdat.Dpar_v = Dpar_smooth_v;
    BSdat.Dperp_v = Dperp_smooth_v;
    BSdat.Dpar = Dpar_smooth;
    BSdat.Dperp = Dperp_smooth;
    BSdat.PDparDperp = PDparDperp_smooth_mean;
    BSdat.PDparDperp_std = PDparDperp_smooth_std;
    BSdat.Diso_v = Diso_smooth_v;
    BSdat.Dratio_v = Dratio_smooth_v;
    BSdat.Diso = Diso_smooth;
    BSdat.Dratio = Dratio_smooth;
    BSdat.PDisoDratio = PDisoDratio_smooth_mean;
    BSdat.PDisoDratio_std = PDisoDratio_smooth_std;
    BSdat.Di_v = Di_smooth_v;
    BSdat.Da_v = Da_smooth_v;
    BSdat.Di = Di_smooth;
    BSdat.Da = Da_smooth;
    BSdat.PDiDa = PDiDa_smooth_mean;
    BSdat.PDiDa_std = PDiDa_smooth_std;

    eval(['save ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/BSdat BSdat'])
    eval(['print -djpeg -loose ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/BSFig'])
end

cd(wd)