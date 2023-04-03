clear all

wd = cd;

DataDir = '/Users/daniel/Dropbox';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
cd(DataDir)


ExpNam = 'AOToct5'; expno = 15; 
signal = 'area';

cd(ExpNam)
for nexp = 1:length(expno)
    eval(['load ' num2str(expno(nexp)) '/FitDat'])

    waterpeak = PP.Npeaks;
    waterpeak = 2;
    nvd = 1;

    PP.peakppm(PP.td1start,waterpeak)
    S = PP.Apeak(:,waterpeak);
    Dpar = FitDat.Dpar(waterpeak);
    Dperp = FitDat.Dperp(waterpeak);
    Diso = (Dpar + 2*Dperp)/3;
    S0 = FitDat.Y0(waterpeak)/AcqDat.Ndir;
    
    ndir = 1:AcqDat.Ndir;
    b_array = reshape(AcqDat.bmat.b,AcqDat.Nb,AcqDat.NDeltab,AcqDat.Ndir);
    S_array = reshape(S,AcqDat.Nb,AcqDat.NDeltab,AcqDat.Ndir,AcqDat.Nvd);
    S_array = squeeze(S_array(:,:,:,nvd));
    Npoints = numel(S_array);
    S_vector = reshape(S_array,Npoints,1);

    NBootStrap = 100;
    ODFsmooth.Pcum = 0;
    rng('default')
    for nBootStrap = 1:NBootStrap
        %DT geometry
        DT.iso = Diso;
        DT.Delta = (Dpar-Dperp)/3/Diso;
        load UniformDistSphereRepulsionN31
        weight = ones(size(theta));
    %     load ZCWfull616
%         load SHREWD_REPfull376
    %    load SHREWD_STEPhemi3321
        ODF.Ndir = length(theta);
        xold = sin(theta).*cos(phi);
        yold = sin(theta).*sin(phi);
        zold = cos(theta);
        steptheta = acos(rand(1,1));
        stepphi = 2*pi*rand(1,1);
        x = xold*cos(steptheta) + zold*sin(steptheta);
        z = zold*cos(steptheta) - xold*sin(steptheta);
        xold = x;
        x = xold*cos(stepphi) - yold*sin(stepphi);
        y = yold*cos(stepphi) + xold*sin(stepphi);                
        ODF.theta = acos(z);
        ODF.phi = atan2(y,x);
        
        figure(1), clf
        plot3(x,y,z,'o')
        axis equal
        xlabel('x')

        %ODF.P = ones(ODF.Ndir,1);
        ODF.P = exp(1*cos(theta)).*weight;
        ODF.P = ODF.P/sum(ODF.P);

        bmat.b = AcqDat.bmat.b;
        bmat.Delta = AcqDat.bmat.Delta;
        bmat.theta = AcqDat.bmat.dir.theta;
        bmat.phi = AcqDat.bmat.dir.phi;
        td1 = numel(bmat.b);

        K.bmat.b = repmat(bmat.b,1,ODF.Ndir);
        K.bmat.Delta = repmat(bmat.Delta,1,ODF.Ndir);
        K.bmat.theta = repmat(bmat.theta,1,ODF.Ndir);
        K.bmat.phi = repmat(bmat.phi,1,ODF.Ndir);
        K.ODF.theta = repmat(ODF.theta',td1,1);
        K.ODF.phi = repmat(ODF.phi',td1,1);

        K.cosbeta = cos(K.bmat.theta).*cos(K.ODF.theta) + sin(K.bmat.theta).*sin(K.ODF.theta).*cos(K.bmat.phi - K.ODF.phi);
        K.P2cosbeta = (3*K.cosbeta.^2 - 1)/2;

        K.Deff = DT.iso.*(1 + 2*K.bmat.Delta.*DT.Delta.*K.P2cosbeta);
        Kernel = exp(-K.bmat.b.*K.Deff);

        %BootStrapPoints = 1:Npoints;
        BootStrapPoints = sort(ceil(Npoints*rand(Npoints,1)));
        Kernel_temp = Kernel(BootStrapPoints,:);
        S_vector_temp = S_vector(BootStrapPoints,1);
        ODF.P = lsqnonneg(Kernel_temp,S_vector_temp/S0);

        Scalc = S0*Kernel_temp*ODF.P;
        Scalc_array = reshape(Scalc,AcqDat.Nb,AcqDat.NDeltab,AcqDat.Ndir);

        figure(2), clf
        plot(1:Npoints,S_vector_temp,'o',1:Npoints,Scalc,'-')

    %     Col = {'r','g','b','k'};
    %     for nDeltab = 1:4
    %         bplot = squeeze(b_array(:,nDeltab,ndir));
    %         Splot = squeeze(S_array(:,nDeltab,ndir));
    %         semilogy(bplot,Splot,['o' Col{nDeltab}])
    %         hold on
    %         Splot = squeeze(Scalc_array(:,nDeltab,ndir));
    %         semilogy(bplot,Splot,['-' Col{nDeltab}])
    %     end
    %     bcalc = max(AcqDat.bmat.b)*linspace(0,1,100);
    %     Sinit = S0*exp(-bcalc*Diso);
    %     semilogy(bcalc,Sinit,'--k')
    %     axis([max(AcqDat.bmat.b)*[-.1 1.1] S0*[1e-2 1.2]])

        %Smooth ODF nodes

        Nphi = 100;
        Ntheta = 100;

        phi = linspace(0,2*pi,Nphi);
        theta = linspace(0,pi,Ntheta);
        dphi = gradient(phi);
        dtheta = gradient(theta);

        [phi,theta] = ndgrid(phi,theta);
        [dphi,dtheta] = ndgrid(dphi,dtheta);
        dOmega = abs(sin(theta).*dtheta.*dphi);

        xnorm = sin(theta).*cos(phi);
        ynorm = sin(theta).*sin(phi);
        znorm = cos(theta);

        ODFsmooth.Nnodes = Nphi*Ntheta;
        ODFsmooth.x = reshape(xnorm,ODFsmooth.Nnodes,1);
        ODFsmooth.y = reshape(ynorm,ODFsmooth.Nnodes,1);
        ODFsmooth.z = reshape(znorm,ODFsmooth.Nnodes,1);

        %Discrete ODF
        ODFdiscrete.x = sin(ODF.theta).*cos(ODF.phi);
        ODFdiscrete.y = sin(ODF.theta).*sin(ODF.phi);
        ODFdiscrete.z = cos(ODF.theta);
        ODFdiscrete.P = ODF.P;

        %Watson distribution smoothing kernel
        ODindex = .02;
        kappa = 1/tan(ODindex*pi/2);

        [K.Xsmooth,K.Xdiscrete] = ndgrid(ODFsmooth.x,ODFdiscrete.x);
        [K.Ysmooth,K.Ydiscrete] = ndgrid(ODFsmooth.y,ODFdiscrete.y);
        [K.Zsmooth,K.Zdiscrete] = ndgrid(ODFsmooth.z,ODFdiscrete.z);

        KernelODFsmooth = exp(kappa*(K.Xsmooth.*K.Xdiscrete + ...
            K.Ysmooth.*K.Ydiscrete + K.Zsmooth.*K.Zdiscrete).^2);

        clear K

        %Smooth ODF amplitude

        ODFsmooth.P = KernelODFsmooth*ODFdiscrete.P;

        ODFsmooth.P = reshape(ODFsmooth.P,[Nphi Ntheta]);
        ODFsmooth.Pcum = ODFsmooth.Pcum + ODFsmooth.P;
        %r = ones(Nphi,Ntheta);
        r = ODFsmooth.P.*ones(Nphi,Ntheta);
        rcum = ODFsmooth.Pcum.*ones(Nphi,Ntheta);
        %r = dOmega.*ones(Nphi,Ntheta);
        Ptheta = sum(ODFsmooth.P.*dOmega,1);
        Pthetacum = sum(ODFsmooth.Pcum.*dOmega,1);

        x = r.*sin(theta).*cos(phi)/max(max(r));
        y = r.*sin(theta).*sin(phi)/max(max(r));
        z = r.*cos(theta)/max(max(r));
        c = ODFsmooth.P;
        c = ones(Nphi,Ntheta,3);
        cnorm = max(max(r));
        cnorm = 1;
        c(:,:,1) = abs(xnorm/cnorm);
        c(:,:,2) = abs(ynorm/cnorm);
        c(:,:,3) = abs(znorm/cnorm);

        xcum = rcum.*sin(theta).*cos(phi)/max(max(rcum));
        ycum = rcum.*sin(theta).*sin(phi)/max(max(rcum));
        zcum = rcum.*cos(theta)/max(max(rcum));

        figure(5), clf
        subplot(1,4,1)
        surf(x,y,z,c)
        shading interp
        axis equal

        subplot(1,4,2)
        surf(xcum,ycum,zcum,c)
        shading interp
        axis equal

        subplot(2,2,2)
        plot(theta(1,:),Ptheta/max(Ptheta),'-',theta(1,:),Pthetacum/max(Pthetacum),'-')

        swh = 2.5e3;
        si = 1024;
        R2 = 10;
        DeltaQ = 550;

        nu = swh*linspace(0,1,si)'; nu = nu - nu(si/2+1);

        [K.nu, K.theta] = ndgrid(nu, theta(1,:));

        K.offset = 2*DeltaQ*(3*cos(K.theta).^2-1)/2;

        Kernel2H = R2.^2./(R2.^2 + 4*(K.offset/2 - K.nu).^2)...
            + R2.^2./(R2.^2 + 4*(-K.offset/2 - K.nu).^2);

        clear K

        I2H = Kernel2H*Ptheta';

        figure(5)
        subplot(2,2,4)
        plot(nu,I2H,'-')
        axis tight off
        pause(.1)
    end

%     Xdat = AcqDat.bmat.b;
%     DirDat = cell(AcqDat.NDeltab*AcqDat.Ndir,1);
%     Dnocalib = zeros(AcqDat.NDeltab*AcqDat.Ndir,PP.Npeaks);
% 
%     PlotInterm = 0; Npeaks = PP.Npeaks;
%     for ndir = 1:(AcqDat.NDeltab*AcqDat.Ndir)
%         fitpoints = (1:AcqDat.Nb) + (ndir-1)*AcqDat.Nb;
%         fitpoints = fitpoints(PP.td1start:length(fitpoints));
%         ExpFit            
%         %GammaFit
%         DirDat{ndir,1} = FitDat;
%         Dnocalib(ndir,:) = FitDat.R(1,:);
%     end 
% %%
%     Dnocalib_array = reshape(Dnocalib(:,1),AcqDat.NDeltab,AcqDat.Ndir);
%     figure(1), clf
%     plot(1:AcqDat.NDeltab,Dnocalib_array,'-o')
%% 
end