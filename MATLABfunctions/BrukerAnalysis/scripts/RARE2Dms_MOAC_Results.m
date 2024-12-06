clear all

wd = cd;

%DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
%DataDir = '/Users/daniel/Dropbox';
%DataDir = '/Users/daniel/Documents/Spaces/Presentations';
DataDir = '/Users/daniel/Dropbox/NMRdata/AVII500/';

%ExpNam = {'AOToct_Eq3'}; expno = 15:10:185;
%ExpNam = {'AOToct_Eq4'}; expno = [15:10:945];
%ExpNam = {'AOToct_Eq5'}; expno = [5:10:485];
%ExpNam = {'AOToct_Eq6'}; expno = [15:10:245]; %expno = 25;
%ExpNam = {'AOToct_Eq7'}; expno = [15:10:245]; %expno = 25;
ExpNam = {'MOACRARE2D_test'}; expno = 60;

fs = 12;
lw = 2;

cd(DataDir)
cd(ExpNam{1})
if exist('expno') == 0
    GetExpnos
end
for nexp = 1:length(expno)
    load([num2str(expno(nexp)) '/NMRacqus'])
        %load([num2str(expno(nexp)) '/PP'])ImagesRaw
        load([num2str(expno(nexp)) '/ImagesRaw'])
        load([num2str(expno(nexp)) '/ImagesMOAC_a'])

        r = ImagesRaw.r;
        [nudim.i,nudim.j,td1] = size(ImagesRaw.S);

        figure(1), clf
        bottom = .7;
        left = .05;
        width = .1;
        height = .1; 
        axh_hist1 = axes('position',[left bottom width height]);
        histdat = sqrt(ImagesMOAC_a.chisq); histdat = reshape(histdat,[numel(histdat),1]);
        hist(histdat)
        axis tight
        hold on
        ylim = get(gca,'YLim'); ylim = abs(diff(ylim))*[-.1 1.1];
        xlim = get(gca,'XLim'); xlim = abs(diff(xlim))*[-.1 .1]+xlim;
        plot([xlim'],[0; 0],'k-')
        set(gca,'XLim',xlim,'YLim',ylim)
        xlabel('rms\chi')
%%        
        ImagesMOAC_BS.S0 = sum(ImagesMOAC_a.w,4);
        ImagesMOAC_BS.MDiso = sum(ImagesMOAC_a.w.*ImagesMOAC_a.iso,4)./ImagesMOAC_BS.S0;
        ImagesMOAC_BS.MR1 = sum(ImagesMOAC_a.w.*ImagesMOAC_a.R1,4)./ImagesMOAC_BS.S0;
        ImagesMOAC_BS.MR2 = sum(ImagesMOAC_a.w.*ImagesMOAC_a.R2,4)./ImagesMOAC_BS.S0;

        ImagesMOAC.S0 = squeeze(mean(ImagesMOAC_BS.S0,5));
        ImagesMOAC.MDiso = squeeze(mean(ImagesMOAC_BS.MDiso,5));
        ImagesMOAC.MR1 = squeeze(mean(ImagesMOAC_BS.MR1,5));
        ImagesMOAC.MR2 = squeeze(mean(ImagesMOAC_BS.MR2,5));
        
        width = .13; height = width*11/8.5;
        bottom = .05; dheight = height+.05;
        left = .05; dleft = width;
        X = r.i*1e3;
        Y = r.j*1e3;
        C = zeros(numel(Y),numel(X),3);

        nslice = 9;
        
        width = 0.4;
        height = 0.4;
        figure(3), clf
        axes('position',[0 .5 width height])
        c.bright = 1*ImagesMOAC.S0(:,:,nslice)'/max(max(max(ImagesMOAC.S0)));
        c.r = ones(size(c.bright)); c.g = ones(size(c.bright)); c.b = ones(size(c.bright));
        C(:,:,1) = c.r.*c.bright; C(:,:,2) = c.g.*c.bright; C(:,:,3) = c.b.*c.bright;
        image(X,Y,C)
        set(gca,'YDir','normal')
        axis equal, axis tight, axis off
        title('S_0')
        
        c.bright = 5*ImagesMOAC.MDiso(:,:,nslice)'/ReconDat.Dmax;
        axes('position',[.5 .5 width height])
        c.r = ones(size(c.bright)); c.g = ones(size(c.bright)); c.b = ones(size(c.bright));
        C(:,:,1) = c.r.*c.bright; C(:,:,2) = c.g.*c.bright; C(:,:,3) = c.b.*c.bright;
        image(X,Y,C)
        set(gca,'YDir','normal')
        axis equal, axis tight, axis off
        title('MD')

        c.bright = ImagesMOAC.MR1(:,:,nslice)'/ReconDat.R1max;
        axes('position',[0 0 width height])
        c.r = ones(size(c.bright)); c.g = ones(size(c.bright)); c.b = ones(size(c.bright));
        C(:,:,1) = c.r.*c.bright; C(:,:,2) = c.g.*c.bright; C(:,:,3) = c.b.*c.bright;
        image(X,Y,C)
        set(gca,'YDir','normal')
        axis equal, axis tight, axis off
        title('R1')

        c.bright = ImagesMOAC.MR2(:,:,nslice)'/ReconDat.R2max;
        axes('position',[.5 0 width height])
        c.r = ones(size(c.bright)); c.g = ones(size(c.bright)); c.b = ones(size(c.bright));
        C(:,:,1) = c.r.*c.bright; C(:,:,2) = c.g.*c.bright; C(:,:,3) = c.b.*c.bright;
        image(X,Y,C)
        set(gca,'YDir','normal')
        axis equal, axis tight, axis off
        title('R2')
return
        figure(4), clf
        axes('position',[0 0 1 1])
        c.bright = 1*ImagesMOAC.S0(:,:,nslice)'/max(max(max(ImagesMOAC.S0)));
        c.r = ones(size(c.bright)); c.g = ones(size(c.bright)); c.b = ones(size(c.bright));
        C(:,:,1) = c.r.*c.bright; C(:,:,2) = c.g.*c.bright; C(:,:,3) = c.b.*c.bright;
        image(X,Y,C)
        set(gca,'YDir','normal')
        axis equal, axis tight, axis off
        title('S_0')

        %Smooth ODF nodes
        load UDSRTriN1000
        ODFsmooth = UDSR;
        TR = TriRep(ODFsmooth.tri, ODFsmooth.x, ODFsmooth.y, ODFsmooth.z);
        Nsubdiv = 1;
        TR=SubdivideSphericalMesh(TR,Nsubdiv);
        ODFsmooth.tri = TR.Triangulation;
        ODFsmooth.x = TR.X(:,1);
        ODFsmooth.y = TR.X(:,2);
        ODFsmooth.z = TR.X(:,3);
        ODFsmooth.N = numel(ODFsmooth.x);
        ODFsmooth.theta = acos(ODFsmooth.z);
        ODFsmooth.phi = atan2(ODFsmooth.y,ODFsmooth.x);

        %Watson distribution smoothing kernel
        ODindex = .05;
        kappa = 1/tan(ODindex*pi/2);

        bottom1 = 0;
        left1 = 0;
        width = 1/nudim.j;
        height = 1/nudim.i;

        for nj = 1:nudim.j    
            for ni = 1:nudim.i
                uDTrec_a.w = squeeze(ImagesMOAC_a.w(ni,nj,nslice,:,:));
                uDTrec_a.iso = squeeze(ImagesMOAC_a.iso(ni,nj,nslice,:,:));
                uDTrec_a.Delta = squeeze(ImagesMOAC_a.Delta(ni,nj,nslice,:,:));
                uDTrec_a.theta = squeeze(ImagesMOAC_a.theta(ni,nj,nslice,:,:));
                uDTrec_a.phi = squeeze(ImagesMOAC_a.phi(ni,nj,nslice,:,:));
                uDTrec_a.R1 = squeeze(ImagesMOAC_a.R1(ni,nj,nslice,:,:));
                uDTrec_a.R2 = squeeze(ImagesMOAC_a.R2(ni,nj,nslice,:,:));

                uDTrec.N = ReconDat.NBS;
                mask = find(uDTrec_a.w>0);
                uDTrec.w = reshape(uDTrec_a.w(mask),[numel(mask) 1])/uDTrec.N;
                uDTrec.iso = reshape(uDTrec_a.iso(mask),[numel(mask) 1]);
                uDTrec.Delta = reshape(uDTrec_a.Delta(mask),[numel(mask) 1]);
                uDTrec.theta = reshape(uDTrec_a.theta(mask),[numel(mask) 1]);
                uDTrec.phi = reshape(uDTrec_a.phi(mask),[numel(mask) 1]);
                uDTrec.R1 = reshape(uDTrec_a.R1(mask),[numel(mask) 1]);
                uDTrec.R2 = reshape(uDTrec_a.R2(mask),[numel(mask) 1]);
                uDTrec.par = uDTrec.iso.*(1 + 2*uDTrec.Delta);
                uDTrec.perp = uDTrec.iso.*(1 - uDTrec.Delta);
                uDTrec.ratio = (1 + 2*uDTrec.Delta)./(1 - uDTrec.Delta);
                uDTrec.N = numel(mask);
                
                [uDTrec.w,indx] = sort(uDTrec.w,'descend');
                uDTrec.iso = uDTrec.iso(indx);
                uDTrec.Delta = uDTrec.Delta(indx);
                uDTrec.theta = uDTrec.theta(indx);
                uDTrec.phi = uDTrec.phi(indx);
                uDTrec.R1 = uDTrec.R1(indx);
                uDTrec.R2 = uDTrec.R2(indx);
                uDTrec.par = uDTrec.par(indx);
                uDTrec.perp = uDTrec.iso(indx);
                uDTrec.ratio = uDTrec.ratio(indx);

                %uDTrec.theta(1:ceil(uDTrec.N*.2)) = NaN;

                mask_ODF = find([uDTrec.ratio>10^1]);
                %mask_ODF = find([uDTrec.ratio>10^.1]);
                %mask_ODF = find([uDTrec.ratio<10^-1]);

                ODFdiscrete.P = uDTrec.w(mask_ODF)/sum(uDTrec.w(mask_ODF));
                ODFdiscrete.x = sin(uDTrec.theta(mask_ODF)).*cos(uDTrec.phi(mask_ODF));
                ODFdiscrete.y = sin(uDTrec.theta(mask_ODF)).*sin(uDTrec.phi(mask_ODF));
                ODFdiscrete.z = cos(uDTrec.theta(mask_ODF));
                

                [K.Xsmooth,K.Xdiscrete] = ndgrid(ODFsmooth.x,ODFdiscrete.x);
                [K.Ysmooth,K.Ydiscrete] = ndgrid(ODFsmooth.y,ODFdiscrete.y);
                [K.Zsmooth,K.Zdiscrete] = ndgrid(ODFsmooth.z,ODFdiscrete.z);

                Kernel = exp(kappa*(K.Xsmooth.*K.Xdiscrete + ...
                    K.Ysmooth.*K.Ydiscrete + K.Zsmooth.*K.Zdiscrete).^2);
                
                indx_powder = find(isnan(uDTrec.theta(mask_ODF)));
                numel(indx_powder)
                if isempty(indx_powder) == 0
                    Kernel(:,indx_powder) = ones(ODFsmooth.N,numel(indx_powder));
                end
                
                Kernel_norm = repmat(sum(Kernel,1),[ODFsmooth.N 1]);
                Kernel = Kernel./Kernel_norm;

                clear K

                %Smooth ODF amplitude
                ODFsmooth.P = Kernel*ODFdiscrete.P;
                ODFsmooth.verts = repmat(ODFsmooth.P,[1 3]).*[sin(ODFsmooth.theta).*cos(ODFsmooth.phi)...
                    sin(ODFsmooth.theta).*sin(ODFsmooth.phi) ...
                    cos(ODFsmooth.theta)];
                ODFsmooth.c = abs([ODFsmooth.x ODFsmooth.y ODFsmooth.z]);


                left = left1 + (ni-1)*width;
                bottom = bottom1 + (nj-1)*height;
                axh_ODFrec = axes('position',[left bottom width height]);
                p = patch('Faces',ODFsmooth.tri,'Vertices',ODFsmooth.verts/max(ODFsmooth.P));
                axis tight, axis square, axis equal
                axis([-1 1 -1 1 -1 1])
                view(0,90)
                xlabel('x'), ylabel('y')
                set(p,'FaceColor','interp','FaceVertexCData',ODFsmooth.c,...
                'EdgeColor','none','LineWidth',1)
                %title('ODF rec')
                axis off
                pause(.1)
            end
        end

        set(gcf, 'PaperPosition', [0 0 20 20],'PaperSize', [20 20]); 
        eval(['print ' num2str(expno(nexp)) '/ODFFig -loose -dpdf'])

 
end
