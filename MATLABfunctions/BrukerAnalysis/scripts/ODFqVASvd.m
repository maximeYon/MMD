clear all

wd = cd;

%DataDir = '/Users/daniel/Dropbox';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
DataDir = '/Users/daniel/Documents/Spaces/Presentations';
cd(DataDir)

ExpNam = 'AOToct5'; expno = 7;
ExpNam = 'RandSampDR1R2'; expno = 65;

td1start = 2;
signal = 'area';
PlotInterm = 0;
NBS = 1000;
minNnodes = 100;
maxNnodes = 150;
%Watson distribution smoothing kernel
ODindex = .05;
kappa = 1/tan(ODindex*pi/2);

cd(ExpNam)
for nexp = 1:length(expno)
    eval(['load ' num2str(expno(nexp)) '/FitDat'])

    waterpeak = PP.Npeaks;
    %waterpeak = 2;

    PP.peakppm(PP.td1start,waterpeak)
    S = PP.Apeak(:,waterpeak);
    Dpar = FitDat.Dpar(waterpeak);
    Dperp = FitDat.Dperp(waterpeak);
    Diso = (Dpar + 2*Dperp)/3;
    DeltaD = (Dpar - Dperp)/3/Diso;
    S0 = FitDat.Y0(waterpeak)/AcqDat.Ndir;
    
    ndir = 1:AcqDat.Ndir;
    b_array = reshape(AcqDat.bmat.b,AcqDat.Nb,AcqDat.NDeltab,AcqDat.Ndir);
    
    for nvd = 1:AcqDat.Nvd
        S_array = reshape(S,AcqDat.Nb,AcqDat.NDeltab,AcqDat.Ndir,AcqDat.Nvd);
        S_array = squeeze(S_array(:,:,:,nvd));

        Npoints = numel(S_array);
        S_vector = reshape(S_array,Npoints,1);
        S0 = max(S_vector);

        load UDSRTriN1000
        ODFsmooth = UDSR;
        ODFsmooth.Pcum = 0;

        if PlotInterm
            figure(2), clf, hold on
        end

        indx_all_array = repmat((1:Npoints)',[1 Npoints]);
        S_array = reshape(S_vector,AcqDat.Nb*AcqDat.NDeltab,AcqDat.Ndir);
        Xin = AcqDat.bmat.b(1:AcqDat.Nb*AcqDat.NDeltab);
        Xin2 = AcqDat.bmat.Delta(1:AcqDat.Nb*AcqDat.NDeltab);
        Pin = [max(S_vector)*1.05 Diso DeltaD]; Funam = 'fDiffVACSY2';
        UB = [max(S_vector)*5 10*Diso 1];
        LB = [max(S_vector)*0.01 .1*Diso -0.5];

        options = optimset('Display','off');
        Pout_array = zeros(NBS,3);

        ODFrec_x = zeros(maxNnodes,NBS);
        ODFrec_y = zeros(maxNnodes,NBS);
        ODFrec_z = zeros(maxNnodes,NBS);
        ODFrec_theta = zeros(maxNnodes,NBS);
        ODFrec_phi = zeros(maxNnodes,NBS);
        ODFrec_P = zeros(maxNnodes,NBS);
        ODFrec_Psmooth = zeros(ODFsmooth.N,NBS);

        ODFrecCell = cell(NBS,6);
        Nnodes = ceil((maxNnodes-minNnodes)*rand(NBS,1) + minNnodes);
        chisq = zeros(NBS,1);

        for nBS = 1:NBS
            eval(['load UDSRTriN' num2str(Nnodes(nBS))])
            ODFrec = UDSR;
            xold = ODFrec.x;
            yold = ODFrec.y;
            zold = ODFrec.z;
            steptheta = acos(2*rand(1,1)-1);
            stepphi = 2*pi*rand(1,1);
            x = xold*cos(steptheta) + zold*sin(steptheta);
            ODFrec.z = zold*cos(steptheta) - xold*sin(steptheta);
            xold = x;
            ODFrec.x = xold*cos(stepphi) - yold*sin(stepphi);
            ODFrec.y = yold*cos(stepphi) + xold*sin(stepphi);                
            ODFrec.theta = acos(ODFrec.z);
            ODFrec.phi = atan2(ODFrec.y,ODFrec.x);

            if PlotInterm
                figure(2)
                trimesh(ODFrec.tri,ODFrec.x,ODFrec.y,ODFrec.z,'EdgeColor','b')
                plot3(1.01*ODFrec.x,1.01*ODFrec.y,1.01*ODFrec.z,'.k')
                axis equal
                xlabel('x')
            end

            ODFrec_x(1:Nnodes(nBS),nBS) = ODFrec.x;
            ODFrec_y(1:Nnodes(nBS),nBS) = ODFrec.y;
            ODFrec_z(1:Nnodes(nBS),nBS) = ODFrec.z;
            ODFrec_theta(1:Nnodes(nBS),nBS) = ODFrec.theta;
            ODFrec_phi(1:Nnodes(nBS),nBS) = ODFrec.phi;

            ODFrecCell{nBS,1} = ODFrec.x;
            ODFrecCell{nBS,2} = ODFrec.y;
            ODFrecCell{nBS,3} = ODFrec.z;
            ODFrecCell{nBS,4} = ODFrec.theta;
            ODFrecCell{nBS,5} = ODFrec.phi;
        end

        tic
        parfor nBS = 1:NBS
        %%    
            %indx_BS = indx_all;
            indx_BS = td1start - 1 + sort(ceil((Npoints-td1start+1)*rand(Npoints,1)));
            %figure(1), clf, plot(indx_all,indx_BS,'-'), return
            indx_BS_array = repmat(indx_BS',[Npoints 1]);

            PAweight = 1e-3 + sum(indx_all_array == indx_BS_array,2);
            %figure(1), clf, plot(indx_all,PAweight,'-'), return
            PAweight_array = reshape(PAweight,AcqDat.Nb*AcqDat.NDeltab,AcqDat.Ndir);

            S_PA = sum(S_array.*PAweight_array,2)./sum(PAweight_array,2);
            %figure(1), clf, plot(1:AcqDat.Nb*AcqDat.NDeltab,S_PA,'o'), return

            Yin = S_PA; 
            %figure(1), clf, semilogy(Xin,Yin,'o'), return    

            Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = Pin; 
            Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,LB./Pnorm,UB./Pnorm,options,Xin2,Pnorm,Xnorm,Ynorm);
            Yout = feval(Funam,Pout,Xin,Xin2); error = Yin - Yout;
            %figure(1), clf, semilogy(Xin,Yin,'o',Xin,Yout,'x'), return
            Pout_array(nBS,:) = Pout;
            S0 = Pout(1); Diso = Pout(2); DeltaD = Pout(3);


            bmat_b = repmat(AcqDat.bmat.b,1,Nnodes(nBS));
            bmat_Delta = repmat(AcqDat.bmat.Delta,1,Nnodes(nBS));
            bmat_theta = repmat(AcqDat.bmat.dir.theta,1,Nnodes(nBS));
            bmat_phi = repmat(AcqDat.bmat.dir.phi,1,Nnodes(nBS));
        %     ODF_theta = repmat(ODFrec_theta(1:Nnodes(nBS),nBS)',AcqDat.td1,1);
        %     ODF_phi = repmat(ODFrec_phi(1:Nnodes(nBS),nBS)',AcqDat.td1,1);
            ODF_theta = repmat(ODFrecCell{nBS,4}',AcqDat.Nb*AcqDat.NDeltab*AcqDat.Ndir,1);
            ODF_phi = repmat(ODFrecCell{nBS,5}',AcqDat.Nb*AcqDat.NDeltab*AcqDat.Ndir,1);

            cosbeta = cos(bmat_theta).*cos(ODF_theta) + sin(bmat_theta).*sin(ODF_theta).*cos(bmat_phi - ODF_phi);
            P2cosbeta = (3*cosbeta.^2 - 1)/2;

            Deff = Diso.*(1 + 2*bmat_Delta.*DeltaD.*P2cosbeta);
            Kernel = exp(-bmat_b.*Deff);

            S_vector_temp = S_vector(indx_BS,1);
            Kernel_temp = Kernel(indx_BS,:);

            ODF_P = lsqnonneg(Kernel_temp,S_vector_temp/S0);
            Scalc = S0*Kernel_temp*ODF_P;
            chisq(nBS,1) = sum((S_vector_temp-Scalc).^2,1)/(AcqDat.Nb*AcqDat.NDeltab*AcqDat.Ndir);

            [Xsmooth,Xrec] = ndgrid(ODFsmooth.x,ODFrecCell{nBS,1});
            [Ysmooth,Yrec] = ndgrid(ODFsmooth.y,ODFrecCell{nBS,2});
            [Zsmooth,Zrec] = ndgrid(ODFsmooth.z,ODFrecCell{nBS,3});

            KernelODFsmooth = exp(kappa*(Xsmooth.*Xrec + ...
                Ysmooth.*Yrec + Zsmooth.*Zrec).^2);

            ODF_Psmooth = KernelODFsmooth*ODF_P;
            ODFrec_Psmooth(:,nBS) = ODF_Psmooth/sum(ODF_Psmooth);

            ODF_Ptemp = zeros(maxNnodes,1);
            ODF_Ptemp(1:Nnodes(nBS),1) = ODF_P;
            ODFrec_P(:,nBS) = ODF_Ptemp;

            if PlotInterm
                Scalc_array = reshape(Scalc,AcqDat.Nb,AcqDat.NDeltab,AcqDat.Ndir);
            end

        %     figure(2), clf
        %     plot(1:Npoints,S_vector_temp,'o',1:Npoints,Scalc,'-')

            nBS

        end
        toc
    %%    
        figure(3), clf
        subplot(2,2,1)
        hist(Pout_array(:,1))
        subplot(2,2,2)
        hist(Pout_array(:,2))
        subplot(2,2,3)
        hist(Pout_array(:,3))
        subplot(2,2,4)
        plot(Pout_array(:,2),Pout_array(:,3),'.')
        %%
        ODFrec_Pav = zeros(ODFsmooth.N,NBS);
        for nBS = 1:NBS
            indx_BS = sort(ceil(NBS*rand(NBS ,1)));
            ODFrec_Pav(:,nBS) = mean(ODFrec_Psmooth(:,indx_BS),2);
        end
        %%
        ODFrec_Pstd = std(ODFrec_Pav,0,2);
        ODFsmooth.Pcum = mean(ODFrec_Psmooth,2);
        ODFrec_Pstd = std(ODFrec_Pav,0,2);
        ODFsmooth.Pcum = mean(ODFrec_Psmooth,2);
        ODFsmooth.Pcum = ODFsmooth.Pcum/sum(ODFsmooth.Pcum)*ODFsmooth.N;

        if PlotInterm
            ODFsmooth.vrec = repmat(ODFsmooth.Prec,[1 3]).*[sin(ODFsmooth.theta).*cos(ODFsmooth.phi)...
            sin(ODFsmooth.theta).*sin(ODFsmooth.phi) ...
            cos(ODFsmooth.theta)];

            ODFsmooth.vreccum = 1/max(ODFsmooth.Pcum)*repmat(ODFsmooth.Pcum,[1 3]).*[sin(ODFsmooth.theta).*cos(ODFsmooth.phi)...
                sin(ODFsmooth.theta).*sin(ODFsmooth.phi) ...
                cos(ODFsmooth.theta)];

            figure(3), clf
            subplot(1,2,1)
            p = patch('Faces',ODFsmooth.tri,'Vertices',ODFsmooth.vrec);
            axis tight, axis square, axis equal
            view(30,30)
            xlabel('x')
            set(p,'FaceColor','interp','FaceVertexCData',ODFsmooth.c,...
            'EdgeColor','k','LineWidth',1)


            subplot(1,2,2)
            p = patch('Faces',ODFsmooth.tri,'Vertices',ODFsmooth.vreccum);
            axis tight, axis square, axis equal
            view(30,30)
            xlabel('x')
            set(p,'FaceColor','interp','FaceVertexCData',ODFsmooth.c,...
            'EdgeColor','k','LineWidth',1)
            pause(.1)
        end

        ODFsmooth.vreccum = repmat(ODFsmooth.Pcum,[1 3]).*[sin(ODFsmooth.theta).*cos(ODFsmooth.phi)...
            sin(ODFsmooth.theta).*sin(ODFsmooth.phi) ...
            cos(ODFsmooth.theta)];
        ODFsmooth.vrecstd = repmat(ODFrec_Pstd,[1 3]).*[sin(ODFsmooth.theta).*cos(ODFsmooth.phi)...
            sin(ODFsmooth.theta).*sin(ODFsmooth.phi) ...
            cos(ODFsmooth.theta)];
        ODFsmooth.c = abs([ODFsmooth.x ODFsmooth.y ODFsmooth.z]);

        figure(6), clf
        axes('position',[0.02 .02 .5*11/8.5 .98],'FontSize',20)
        p = patch('Faces',ODFsmooth.tri,'Vertices',ODFsmooth.vreccum);
        axis tight, axis square, axis equal
        axis(2.5*[-1 1 -1 1 -1 1])
        view(30,30)
        set(p,'FaceColor','interp','FaceVertexCData',ODFsmooth.c,...
        'EdgeColor','k','LineWidth',1)
        set(gca,'LineWidth',5,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[])
        %axis off
        eval(['print ' DataDir '/' ExpNam '/' num2str(expno) '/ODF' num2str(nvd) ' -djpeg -r300'])

        subplot(2,3,6), cla
        p = patch('Faces',ODFsmooth.tri,'Vertices',ODFsmooth.vrecstd);
        axis tight, axis square, axis equal
        axis([xlim xlim xlim])
        view(30,30)
        xlabel('x')
        set(p,'FaceColor','interp','FaceVertexCData',ODFsmooth.c,...
        'EdgeColor','k','LineWidth',1)

        figure(4), clf
        plot(Nnodes,sqrt(chisq),'o')
        hold on
        %plot([min(Nnodes); max(Nnodes)],1/SNR*[1; 1],'k-')
    end

    %% 
end