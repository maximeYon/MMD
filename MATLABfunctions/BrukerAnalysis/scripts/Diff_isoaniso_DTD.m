clear all

wd = cd;

%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/Dropbox';
DataDir = '/Users/daniel/Dropbox/NMRdata/AVII500';

ExpNam = 'Yeast_tripgse'; expno = [52:58 67 69 73 75 77]; expno = [73 75 77]; expno = 73;

td1start = 1;
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
    
    b = AcqDat.bmat.b;
    td1 = AcqDat.td1;
    bt = AcqDat.bmat;
    signal = S;
    ind = td1start:td1;
        
    opt.dtd.do_plot = 1;
    opt = mdm_opt(opt);
    opt = dtd_opt(opt);
    
    xps.b = b;
    xps.n = td1;
    xps.bt = [bt.xx bt.yy bt.zz sqrt(2)*[bt.xy bt.xz bt.yz]];

    [a_ind,~] = ndgrid(1:AcqDat.Nb,1:AcqDat.Ndir);
    xps.a_ind = reshape(a_ind,xps.n,1); % Same a_ind means same b-tensor size and shape

    xps_old = xps;
    xps = mdm_xps_calc_btpars(xps_old);
    
%     figure(1), clf
%     plot3(xps.b_l,xps.b_s,abs(xps.b_lambdazzvec(:,3)),'.')
%     return

    %dtd = [4 2e-9 2e-9 0 0 .4 2e-11 2e-11 0 0 .1 2e-9 1e-11 0 0 .25 2e-9 1e-11 pi/4 0 .25]';
    Dpar = 2e-9; Dperp = 1e-10; theta = 45/180*pi;
    dtd_cross = [3 Dpar Dperp theta 0 1/3 Dpar Dperp theta 120/180*pi 1/3 Dpar Dperp theta 240/180*pi 1/3]';
    [dtd_nx6,w] = dtd_dist2nx6w(dtd_cross);
    s0 = ones(1,numel(w))*w;
    dt1x6 = (dtd_nx6'*w)'/s0;

    Dpar_fat = dt1x6(3)
    Dperp_fat = dt1x6(1)
    dtd_fat = [1 Dpar_fat Dperp_fat 0/180*pi 0/180*pi 1]';
    
    Dpar_bi = .8e-9;
    Dperp_bi = 1e-10;
    f_s = 1 - (Dpar_fat-Dperp_fat)/(Dpar_bi-Dperp_bi)
    D_s = Dpar_bi - (Dpar_bi-Dpar_fat)
    dtd_bi = [2 Dpar_bi Dperp_bi 0/180*pi 0/180*pi 1-f_s  D_s D_s 0/180*pi 0/180*pi f_s]';
     
    dtd = dtd_bi;
    m = dtd_dtd2m(dtd,opt);
    signal = dtd_1d_fit2data(m, xps);

    [n,par,perp,theta,phi,w] = dtd_dist2par(dtd);
                                
    if n > 0
        s0 = ones(1,n)*w;
        iso_v = (par + 2*perp)/3;
        miso = iso_v'*w/s0;
        viso = (iso_v-miso)'.^2*w/s0;
        aniso_v = par - perp;
        maniso = aniso_v'*w/s0;
        vaniso = (aniso_v-maniso)'.^2*w/s0;
        vlambda_v = 2/9*aniso_v.^2;
        mvlambda = vlambda_v'*w/s0;
        vvlambda = (vlambda_v-mvlambda)'.^2*w/s0;

        [dtd_nx6,w] = dtd_dist2nx6w(dtd);
        dt1x6 = (dtd_nx6'*w)'/s0
        dt3x3 = tm_1x6_to_3x3(dt1x6);

        dt = tm_t2tpars(dt3x3);
    end
    
mu2iso = viso;
mu2aniso = 2/5*mvlambda;
mu2tot = mu2iso + mu2aniso;
ciso = mu2iso./miso.^2;
caniso = mu2aniso./miso.^2;
   
    
    id = xps.a_ind;
    [~,~,id_ind] = unique(id, 'rows');

    % get rid of NaNs
    tmp = sum(isnan(id),2) > 0;
    c_list = unique(id_ind(~tmp)); 
    n = numel(c_list);

    signal_pa = zeros(n,1);

    for c = c_list'
        signal_pa(c == c_list,1) = nanmean(signal(id_ind == c,1),1);
    end
    
    % Average fields in xps
    clear xps_pa
    xps = rmfield(xps, 'n');
    f = fieldnames(xps);
    for i = 1:numel(f)
        for c = c_list'

            % allow text fields to be just copied
            if (all(ischar(xps.(f{i}))))
                xps_pa.(f{i}) = xps.(f{i});
                continue; 
            end

            try
                xps_pa.(f{i})(c == c_list,:) = mean(xps.(f{i})(id_ind == c, :), 1);
            catch
                error('failed powder averaging field %s', f{i});
            end
        end
    end
    xps_pa.n = n;

    figure(1), clf
    plot3(xps_pa.b_s,xps_pa.b_l,log10(signal_pa),'x')
    view(150,10)
    set(gca,'ZLim',[-2 1.1])
    return

%     m = dtd_1d_data2fit(signal, xps, opt, ind);
    %%
    dmin = opt.dtd.dmin;
    dmax = opt.dtd.dmax;
    ratiomax = dmax/dmin;

    nx = 30; ny = nx;
    sigma_x = .1;
    sigma_y = sigma_x;
    xmin = log10(dmin)-3*sigma_x;
    xmax = log10(dmax)+3*sigma_x;
    ymin = log10(1/ratiomax)-1*sigma_y;
    ymax = log10(ratiomax)+1*sigma_y;

    x = linspace(xmin,xmax,nx)';
    dx = x(2)-x(1);
    y = linspace(ymin,ymax,ny)';
    dy = y(2)-y(1);
    [xx,yy] = ndgrid(x,y);

    np = nx*ny;
    
    dtd = dtd_m2dtd(m);
    [n,par,perp,theta,phi,w] = dtd_dist2par(dtd);
    if n>0
        logiso = log10((par + 2*perp)/3);
        logratio = log10(par./perp);

        x0 = logiso(:);
        y0 = logratio(:);
        nx0 = n;
        xx_k = repmat(xx(:),[1 nx0]);
        yy_k = repmat(yy(:),[1 nx0]);
        x0_k = repmat(x0',[nx*ny 1]);
        y0_k = repmat(y0',[nx*ny 1]);
        k = 1/(sigma_x*sqrt(2*pi)).*exp(-(xx_k-x0_k).^2/(2*sigma_x^2)).*...
            1/(sigma_y*sqrt(2*pi)).*exp(-(yy_k-y0_k).^2/(2*sigma_y^2));

        p = k*w;   
    end
                
    z = reshape(p,[nx ny]);

    zprojx = sum(z,2);
    zprojy = sum(z,1);

    nclevels = 12;
    clevels = max(z(:))*linspace(0,1,nclevels+2);
    clevels = clevels(2:(nclevels+1));

    figsize = 2*3.3*[1.618 1];
    figaspect = figsize(1)/figsize(2);

    fs = 2*7;
    lw = 2*1;

    figure(2), clf
    set(gcf, 'PaperUnits','inches', 'PaperPosition', 1*[0 0 figsize],'PaperSize', figsize);

    left = .7;
    bottom = .6;
    height = .3;
    width = height/figaspect;
    proj_height = height*.3;
    proj_width = proj_height/figaspect;

    axes('position',[left bottom width height])
    contour(x,y,z',clevels,'k','LineWidth',.25*lw)
    set(gca,'XLim',[xmin xmax], 'YLim',[ymin ymax],'YAxisLocation','right',...
    'XTick',[-11:-8],'YTick',-2:2,'TickDir','out','TickLength',.03*[1 1],...
    'FontSize',fs,'LineWidth',lw)
    axis square
    xlabel('log(D_{iso} / m^2s^-^1)','FontSize',fs)
    ylabel('log(D_{||} / D_{\perp})','FontSize',fs)

    axes('position',[left bottom+height width proj_height])
    plot(x,zprojx/max(zprojx),'k-','LineWidth',lw)
    set(gca,'XLim',[xmin xmax], 'YLim',[-.3 1.1])
    axis off

    axes('position',[left-proj_width bottom proj_width height])
    plot(zprojy/max(zprojy),y,'k-','LineWidth',lw)
    set(gca,'YLim',[ymin ymax], 'XLim',[-.3 1.1],'XDir','reverse')
    axis off

                    return
    
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