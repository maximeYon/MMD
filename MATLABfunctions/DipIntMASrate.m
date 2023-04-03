clear all

cputime1=cputime;

td = 1*1024; %no of points in time domain signal
swh = 100e3; %spectral width [Hz]

offset = 0;
Deltanudip = 10e3;
omegaRvector = [0 .125 .25 .5 1]*Deltanudip*2*pi;
%omegaRvector = linspace(0,1,500)*Deltanudip*2*pi;



for nomegaR = 1:length(omegaRvector)
    omegaR = omegaRvector(nomegaR);


    dw = 1/swh; %dwell time [s]
    si = 4*td; %no of points in spectrum
    sr = swh/si; %digital resolution in spectrum [Hz]
    aq = td*dw; %acquisition time [s]
    T2 = aq/5; %decay time of the time domain signal, T2 [s]

    t = aq*linspace(0,1,td); %time vector

    nu = swh*linspace(0,1,si); %frequency vector
    nu = nu - nu(si/2+1); %set the proper zero frequency

    Ntimes = td;

    dt = dw;

    theta = linspace(0,-pi/2,200);
    phi = zeros(size(theta));
    weight = abs(sin(theta));
    weight = weight/sum(weight);

    %load LEBhemi151
    %load LEBhemi385
    %load STEPhemi900
    %load SHREWD_STEPhemi66
    %load SHREWD_STEPhemi378
    %load SHREWD_STEPhemi946
    %load SHREWD_STEPhemi3321
    %load REPfull376
    %load SHREWD_REPfull376
    %load ZCWfull377
    load ZCWfull616


    alpha = acos(sqrt(1/3)) - 0*pi/180;

    phi0 = phi;
    theta0 = theta;
    Nphi = length(phi0);

    z0 = cos(theta0);
    x0 = sin(theta0).*cos(phi0);
    y0 = sin(theta0).*sin(phi0);

    z0prim = cos(theta0);
    x0prim = sin(theta0).*cos(phi0);
    y0prim = sin(theta0).*sin(phi0);

    if omegaR ~= 0
        x0 = x0prim*cos(alpha) - z0prim*sin(alpha);
        y0 = y0prim;
        z0 = z0prim*cos(alpha) + x0prim*sin(alpha);
    end


    [x0Array,tArray] = ndgrid(x0,t);
    [y0Array,tArray] = ndgrid(y0,t);
    [z0Array,tArray] = ndgrid(z0,t);
    [weightArray,tArray] = ndgrid(weight,t);

    gammaArray = omegaR*tArray;

    xprim = x0Array*cos(alpha) + z0Array*sin(alpha);
    yprim = y0Array;
    zprim = z0Array*cos(alpha) - x0Array*sin(alpha);

    clear x0Array y0Array z0Array

    xyprim = xprim + i*yprim;
    xyprim = xyprim.*exp(i*gammaArray);
    xprim = real(xyprim);
    yprim = imag(xyprim);
    zprim = zprim.*ones(size(gammaArray));

    x = xprim*cos(alpha) - zprim*sin(alpha);
    y = yprim;
    z = zprim*cos(alpha) + xprim*sin(alpha);

    clear xprim yprim zprim xyprim

    %rprim = sqrt(xprim.^2 + yprim.^2 + zprim.^2);
    %r = sqrt(x.^2 + y.^2 + z.^2);

    thetaArray = acos(z);

    thetaSphere = linspace(0,pi,100);
    phiSphere = linspace(0,2*pi,200);

    [thetaSphere,phiSphere] = ndgrid(thetaSphere,phiSphere);

    P2cosSphere = 0.5*(3*cos(thetaSphere).^2-1);

    xSphere = .99*sin(thetaSphere).*cos(phiSphere);
    ySphere = .99*sin(thetaSphere).*sin(phiSphere);
    zSphere = .99*cos(thetaSphere);

    figure(1), clf
    surf(xSphere,ySphere,zSphere,P2cosSphere)
    shading('interp'), axis('square')
    axis off
    hold on
    plot3(2*[-1 1],[0 0],[0 0],'-k','LineWidth',3)
    plot3([0 0],2*[-1 1],[0 0],'-k','LineWidth',3)
    plot3([0 0],[0 0],2*[-1 1],'-k','LineWidth',3)
    axis([-1 1 -1 1 -1 1])
    view(-150,30)
    caxis([-1 1])

    plot3(x0',y0',z0','ok')
    plot3(x',y',z','-k')
    hold off

    colormap('default')
    colormap('hot')
    hot = colormap;
    hot = hot(20:64,:);
    cool = [hot(:,3) hot(:,2) hot(:,1)];
    hotcool = [cool; flipud(hot)];
    colormap(hotcool)
    colorbar
    colorbar('YTick',[-.5 0 .5])

    clear x y z

    P2costheta = 0.5*(3*cos(thetaArray).^2-1);


    omegaArray = 2*pi*Deltanudip*P2costheta;
    omegaArray = [zeros(Nphi,1) omegaArray(:,1:Ntimes-1)];
    phaseArray = cumsum(omegaArray*dt,2);

%     figure(2), clf
%     subplot(2,1,1)
%     plot(omegaR*t/2/pi,omegaArray/2/pi,'k-','LineWidth',.5)
%     axis([0 max(omegaR*t/2/pi) 1.5*Deltanudip*[-.75 1]])
%     ylabel('\Delta\nu_d / Hz','FontSize',20)
%     %legend('\nu_A','\nu_A')
%     set(gca,'LineWidth',2)
% 
%     subplot(2,1,2)
%     plot(omegaR*t/2/pi,phaseArray,'k-','LineWidth',.5)
%     %hold on
%     % plot(t,-phaseArray,'--','LineWidth',3)
%     % hold off
%     axis([0 max(omegaR*t/2/pi) 1.5*max(max(phaseArray))*[-1 1]])
%     ylabel('\theta / rad','FontSize',20)
%     set(gca,'LineWidth',2)
%     xlabel('\omega_R\itt\rm/2\pi','FontSize',20)
%     %print DipIntMASplot -djpeg -r300

    SignalArray = exp(-i*phaseArray - i*2*pi*offset.*tArray) + exp(i*phaseArray - i*2*pi*offset.*tArray);
    %return
    SArray = weightArray.*SignalArray.*exp(-tArray/T2);

    clear omegaArray phaseArray SignalArray weightArray tArray

    S = sum(SArray(:,:),1);

    clear SArray

    S(:,1) = .5*S(:,1); %multiply first point by .5 to get zero baseline of spectrum
    I = ifft(S,si,2); %fast Fourier transform
    I = ifftshift(I,2); %swap left and right parts of spectrum


    figure(3), clf
    subplot(2,1,1)
    plot(t,real(S),'-',t,imag(S),'-')
    xlabel('time / s')
    legend('real','imag')
    %title(['offsetA=' num2str(offsetA) 'Hz offsetB=' num2str(offsetB) 'Hz k=' num2str(k) 'Hz T_2=' num2str(T2) 's'])

    subplot(2,1,2)
    plot(nu,real(I),'-')
    xlabel('frequency / Hz')
    axis( [offset + 3*Deltanudip*[-1 1] max(real(I))*[-.1 1.1]])
    set(gca,'XDir','reverse')
    title(['swh=' num2str(swh) 'Hz sr=' num2str(sr,3) 'Hz si=' num2str(si) ' aq=' num2str(aq) 's dw=' num2str(dw,3) 's td=' num2str(td)])

    % figure(3), clf
    % subplot(3,1,1)
    % plot(t,thetaArray)
    % ylabel('\theta')
    % title(['\alpha=' num2str(alpha)])
    % subplot(3,1,2)
    % plot(t,omegaArray/2/pi)
    % ylabel('\nu')
    % subplot(3,1,3)
    % plot(t,phaseArray)
    % xlabel('t'), ylabel('\theta')

    Iall(nomegaR,:) = real(I);
    pause(.1)
end

Iallplot = zeros(size(Iall));
Imax = max(max(Iall));
for nomegaR = 1:length(omegaRvector)
    Iallplot(nomegaR,:) = .9*Iall(nomegaR,:)/max(Iall(nomegaR,:)) + 1*(nomegaR-1);
end

nuplot = nu/Deltanudip;

figure(3), clf
axes('position',[0 .2 1 .8])
plot(nuplot,Iallplot,'-k','LineWidth',3)
xlabel('(\it\nu\rm -  \it\nu\rm_A) / \itD\rm_A_X','FontSize',25)
axis([4*[-1 1] [-.2 length(omegaRvector)+.2]])
set(gca,'XDir','reverse','FontSize',20,'LineWidth',2,'YTick',[],'XTick',[-2:1:2],'box','off','TickDir','out')

for nomegaR = 1:length(omegaRvector)
    text(-1.2,nomegaR-.75,num2str(omegaRvector(nomegaR)/(Deltanudip*2*pi),2),'FontSize',20)
end
print DipIntMASrate -deps

return
Iallplot = zeros(size(Iall));
for nomegaR = 1:length(omegaRvector)
    Iallplot(nomegaR,:) = log10(Iall(nomegaR,:)/max(Iall(nomegaR,:)));
end

logthresh = -1;
Iallplot(find(logthresh>Iallplot)) = logthresh;

nuRplot = omegaRvector/2/pi/Deltanudip;

figure(4), clf
surf(nuplot,nuRplot,Iallplot)
shading('interp')
view(0,90)
colormap('hot')
axis('tight')
set(gca,'XDir','reverse','FontSize',20,'LineWidth',2,'YTick',[0 .5 1],'Xlim',1.5*[-1 1],'box','off','TickDir','out')
xlabel('(\it\nu\rm -  \it\nu\rm_A) / \itD\rm_A_X','FontSize',25)
ylabel('\it\nu\rm_R / \itD\rm_A_X','FontSize',25)

%print DipIntMASratecol -depsc

cputime2=cputime-cputime1
