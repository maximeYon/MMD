clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/opt/topspin2/data/DT/nmr';
DataDir = '/Users/daniel/Dropbox';

ExpNam = 'Yeast_tripgse'; expno = [52:58 67 69 73 75 77]; expno = [73 75 77];


for nexp = 1:length(expno)
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/NMRacqus'])
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/FitDat'])
    eval(['load ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/BSdat'])
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
    Yin = S_PA/S0;
    %figure(1), clf, plot(1:AcqDat.Nb,Yin,'o'), return
    
    Xin = AcqDat.bmat.b(1:AcqDat.Nb);
    Xin2 = AcqDat.bmat.Delta(1:AcqDat.Nb);

    Dpar_v = BSdat.Dpar_v;
    Dperp_v = BSdat.Dperp_v;
            
    Diso_v = (Dpar_v + 2*Dperp_v)/3;
    DDelta_v = (Dpar_v - Dperp_v)/3./Diso_v;
    %figure(1), clf, semilogx(Diso_v,DeltaD_v,'o'), return

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
    %figure(1), clf, plot((1:AcqDat.Nb)',Kernel,'-'), set(gca,'YLim',[-.1 1.1]), return
    %figure(1), clf, plot(K_bDDelDel,Kernel,'o'), return

    Ycalc = Kernel*BSdat.PDparDperp;
    %figure(1), clf, plot((1:AcqDat.Nb)',Yin,'o',(1:AcqDat.Nb)',Ycalc,'-'), set(gca,'YLim',[-.1 1.1]), return
%%    
%     indx = find(BSdat.PDparDperp>1e-1*max(BSdat.PDparDperp));
%     Dpar = BSdat.Dpar_v(indx);
%     Dperp = BSdat.Dperp_v(indx);
%     Diso = Diso_v(indx);
%     PD = BSdat.PDparDperp(indx);
%     ND = 1024;
%     D = max([Dpar; Dperp])*linspace(-.05,1.05,ND);
%     [K_D,K_Dpar] = ndgrid(D,Dpar);
%     [K_D,K_Dperp] = ndgrid(D,Dperp);
%     Kernel = 1/2./sqrt((K_D-K_Dperp).*(K_Dpar-K_Dperp));
%     indx = find(K_D > max(cat(3,K_Dperp,K_Dpar),[],3));
%     Kernel(indx) = 0;
%     indx = find(K_D < min(cat(3,K_Dperp,K_Dpar),[],3));
%     Kernel(indx) = 0;
%     Dfreq = linspace(0,1,ND)';
%     Dfreq = Dfreq - Dfreq(ND/2+1);
%     PDft = fftshift(fft(Kernel,ND,1),1);
%     lb = 5;
%     lbfun = exp(-(lb*Dfreq*pi).^2);
%     lbfunarray = repmat(lbfun,1,length(PD));
%     PDft = lbfunarray.*PDft;
%     Kernel = real(ifft(ifftshift(PDft,1),ND,1));
%     figure(1), clf, plot(D,Kernel,'-'), return
%%

    fs = 12;
    lw = 2;

    figure(1), clf
    axes('position',[-.02 .1 .9 .2],'FontSize',fs*.8)
    plot(PP.ppm,PP.Ispec,'k','LineWidth',lw)
    set(gca,'XDir','reverse','YTick',[],'LineWidth',1.5*lw,'Box','off','TickDir','out','Ycolor',[1 1 1])
    axis('tight')
    ylim = get(gca,'YLim');
    ylim = .05*diff(ylim)*[-1 1] + ylim;
    set(gca,'YLim',ylim)
    xlabel(['\delta(' NMRacqus.nuc1 ') / ppm'],'FontSize',fs)
    hold on
    plot(PP.peakppm(PP.td1start,:),PP.Ipeak(PP.td1start,:),'bo','LineWidth',lw)

    [X,Y] = fSchmidt(AcqDat.bmat.dir.x,AcqDat.bmat.dir.y,AcqDat.bmat.dir.z);
    latitude.theta = pi/180*[30:30:150 179];
    latitude.phi = linspace(0,2*pi,100);
    [latitude.phi,latitude.theta] = ndgrid(latitude.phi,latitude.theta);
    latitude.z = cos(latitude.theta);
    latitude.x = sin(latitude.theta).*cos(latitude.phi);
    latitude.y = sin(latitude.theta).*sin(latitude.phi);
    [latitude.X,latitude.Y] = fSchmidt(latitude.x,latitude.y,latitude.z);
    longitude.theta = pi/180*linspace(30,180,100);
    longitude.phi = pi/180*[30:30:360];
    [longitude.theta,longitude.phi] = ndgrid(longitude.theta,longitude.phi);
    longitude.z = cos(longitude.theta);
    longitude.x = sin(longitude.theta).*cos(longitude.phi);
    longitude.y = sin(longitude.theta).*sin(longitude.phi);
    [longitude.X,longitude.Y] = fSchmidt(longitude.x,longitude.y,longitude.z);

    axes('position',[.79 .12 .18 .18],'FontSize',fs*1)
    plot(X,Y,'ko','MarkerFaceColor',[0 0 0])
    hold on
    plot(latitude.X,latitude.Y,'b-')
    plot(longitude.X,longitude.Y,'b-')
    axis tight equal off
    title(['Ndir=' num2str(AcqDat.Ndir)])

    bottom = .45;
    left = .15;
    width = .3;
    Imin = 2e-2;

    Xval = 1e-9*reshape(AcqDat.bmat.iso(1:AcqDat.Nb),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    Yval = 1e-9*reshape(AcqDat.bmat.aniso(1:AcqDat.Nb),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    axh1 = axes('position',[left bottom width .5]);
    Zval = reshape(log10(Ycalc),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    Zval(Zval<log10(Imin)) = NaN;
    ph2 = plot3(Xval,Yval,Zval,'k-');
    hold on
    ph3 = plot3(Xval',Yval',Zval','k-');
    Zval = reshape(log10(Yin),sqrt(AcqDat.Nb),sqrt(AcqDat.Nb));
    Zval(Zval<log10(Imin)) = NaN;
    mh1 = plot3(Xval,Yval,Zval,'bo');
    view(140,10)

    axis([max(AcqDat.bmat.iso)*1e-9*[-.1 1.1 -.1 1.1] log10(Imin) .1])
    set([mh1],'LineWidth',.1*lw,'MarkerSize',3*lw,'MarkerFaceColor',[0 0 1])
    set([ph2 ph3],'LineWidth',1*lw)
    Itick = [.01 .1 1]; ItickLabel = {'0.01', '0.1', '1'};
    set(axh1,'ZTick',log10(Itick),'ZTickLabel',ItickLabel)
    set([axh1 ],'LineWidth',lw,'FontSize',fs,'TickDir','out','Box','off')
    xlabel('b_{iso} / 10^9 sm^-^2'), ylabel('b_{aniso} / 10^9 sm^-^2'), zlabel('S/S_0')
    title(['Nb=' num2str(AcqDat.Nb)])

    left = .60;
    
    Xval = log10(BSdat.Diso);
    Yval = log10(BSdat.Dratio);
    Zval = reshape(BSdat.PDisoDratio,numel(Xval),numel(Yval))';
%     Xval = log10(BSdat.Di);
%     Yval = log10(BSdat.Da);
%     Zval = reshape(BSdat.PDiDa,numel(Xval),numel(Yval))';
    ZprojX = sum(Zval,1);
    ZprojY = sum(Zval,2);
    axes('position',[left bottom width width*11/8.5])
%     imagesc(Xval,Yval,Zval)
%     axis tight, shading flat, colormap('hot')
    minclevel = .05;
    maxclevel = .9;
    Nclevel = 10;
    clevels = max(max(Zval))*linspace(minclevel,maxclevel,Nclevel);
    %clevels = max(max(Zval))*logspace(log10(minclevel),log10(maxclevel),Nclevel);
    contour(Xval,Yval,Zval,clevels,'k','LineWidth',.5*lw);
    set(gca,'LineWidth',lw,'FontSize',fs,'YDir','normal','YAxisLocation','right','TickDir','out','TickLength',.02*[1 1])
    xlabel('log(D_{iso} / m^2s^-^1)'), ylabel('log(D_{par}/D_{perp})')
    axes('position',[left bottom+width*11/8.5+.02 width .12]);
    plot(Xval,ZprojX,'k-','LineWidth',lw)
    axis tight off
    ylim = get(gca,'YLim'); ylim = abs(ylim(2)-ylim(1))*[-.05 1.05]; set(gca,'YLim',ylim)
    axes('position',[left-.11 bottom .12*8.5/11 width*11/8.5]);
    plot(ZprojY,Yval,'k-','LineWidth',lw)
    set(gca,'XDir','reverse')
    axis tight off
    xlim = get(gca,'XLim'); xlim = abs(xlim(2)-xlim(1))*[-.05 1.05]; set(gca,'XLim',xlim)
    
    

    %
%    eval(['save ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/BSdat BSdat'])
%    eval(['print -djpeg -loose ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/BSFig'])
    delete([DataDir '/' ExpNam '/' num2str(expno(nexp)) '/ReportFig.*'])
    eval(['print -depsc -loose ' DataDir '/' ExpNam '/' num2str(expno(nexp)) '/ReportFig'])
end

cd(wd)