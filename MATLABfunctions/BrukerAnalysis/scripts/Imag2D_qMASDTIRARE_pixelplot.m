clear allwd = cd;%DataDir = '/opt/topspin2/data/DT/nmr';DataDir = '/Users/daniel/NMRdata/AVII500/DT';%ExpNam = {'DTItest'}; expno = 125; %expno = 347:360;%ExpNam = {'NafionDTI'}; expno = 31; %expno = [26:28 31:38]; %ExpNam = {'DTI_HHnerve'}; %expno = 415;%ExpNam = {'MIC5_10mm_gradcalib'}; expno = 130;%ExpNam = {'isoaniso_test'}; expno = 104;%ExpNam = {'qMASdtirare2d_test'}; expno = 10;ExpNam = {'uFA_RARE'}; expno = [113:-1:101];%expno = 112; Pixel.x = [67; 67; 74; 67]; Pixel.y = [20; 40; 54; 99]; lb = 4e-4;%expno = 102; Pixel.x = [67; 70; 67; 93]; Pixel.y = [37; 76; 20; 86]; lb = 2e-4;%expno = 111; Pixel.x = [67; 70; 67; 93]; Pixel.y = [37; 76; 20; 86]; lb = 2e-4;ExpNam = {'SjolundOpt'}; expno = 13;Pixel.x = [36 36 36 45]; Pixel.y = [16 35 45 50]; lb = 150e-6; si = 64; si1 = si;fs = 4*11/3.33; lw = .5*11/3.33; tl = .03;%si = 128; si1 = si;td1start = 2;%td1end = 16;thresh = .1;MakeFit = 0;fontsize = 15;cd(DataDir)%GetExpnamsfor ndir = 1:length(ExpNam)    ExpNam{ndir}    cd(ExpNam{ndir})    if exist('expno') == 0        GetExpnos    end    for nexp = 1:length(expno)        expno(nexp)        ConvertAcqus = 'Y';    ConvertProcs = 'N';    MakeTextfile = 'N';        if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0            ConvertAcqus = 'Y';        end        res = fExpnoInfo2(ConvertAcqus,MakeTextfile,ConvertProcs,expno(nexp));        eval(['load ' num2str(expno(nexp)) '/NMRacqus'])        if any(strcmp(NMRacqus.pulprog,{'DT_qMASdtirare2d'})) == 1                        RARE2Dtd1                        %figure(1), clf, imagesc(abs(Itd1(:,:,2))'), axis square, return                        Imax = abs(Itd1);            Imax(nudim.i/2+(0:2),nudim.j/2+(0:2),2:td1) = 0;            Imax = max(reshape(Imax,numel(Imax),1));            fid = fopen([num2str(expno(nexp)) '/rx.txt']);            ramp.x = fscanf(fid,'%f');            fclose(fid);            fid = fopen([num2str(expno(nexp)) '/ry.txt']);            ramp.y = fscanf(fid,'%f');            fclose(fid);            fid = fopen([num2str(expno(nexp)) '/rz.txt']);            ramp.z = fscanf(fid,'%f');            fclose(fid);            fid = fopen([num2str(expno(nexp)) '/riso.txt']);            ramp.iso = fscanf(fid,'%f');            fclose(fid);            fid = fopen([num2str(expno(nexp)) '/qMASx.txt']);            Gmod.x = fscanf(fid,'%f');            fclose(fid);            fid = fopen([num2str(expno(nexp)) '/qMASy.txt']);            Gmod.y = fscanf(fid,'%f');            fclose(fid);            fid = fopen([num2str(expno(nexp)) '/qMASz.txt']);            Gmod.z = fscanf(fid,'%f');            fclose(fid);            fid = fopen([num2str(expno(nexp)) '/qMASr.txt']);            Gmod.r = fscanf(fid,'%f');            fclose(fid);            G.x = ramp.x*Gmax.*NMRacqus.cnst1/100;            G.y = ramp.y*Gmax.*NMRacqus.cnst2/100;            G.z = ramp.z*Gmax.*NMRacqus.cnst3/100;            G.r = sqrt(G.x.^2 + G.y.^2 + G.z.^2);            G.iso = ramp.iso*Gmax.*NMRacqus.cnst3/100;            G.max = max([G.iso; G.r]);            Gnorm.x = G.x./G.r;            Gnorm.y = G.y./G.r;            Gnorm.z = G.z./G.r;                        Ndt = length(Gmod.x);            tau = NMRacqus.d3;            dt = tau/Ndt;            F.r = cumsum(Gmod.r*dt);            b.r = sum(F.r.^2*dt);            F.x = cumsum(Gmod.x*dt);            F.y = cumsum(Gmod.y*dt);            F.z = cumsum(Gmod.z*dt);            F.iso = sqrt(F.x.^2 + F.y.^2 + F.z.^2);            b.iso = sum(F.iso.^2*dt);            b = 2*b.iso*gamma^2*(G.r + G.iso).^2;%             indx_iso = find(isnan(Gnorm.x) == 1);%             Ndir = indx_iso(2) - indx_iso(1);            tol = 1e-3;            minb = min(b);            indx_b0 = find([b>minb*(1-tol) & b<minb*(1+tol)]);            Ndir = length(indx_b0);            figure(1), clf            subplot(2,2,1)            plot(1:td1,b)            subplot(2,2,2)            plot(1:length(Gmod.x),Gmod.x,'r-',1:length(Gmod.x),Gmod.y,'g-',...                1:length(Gmod.x),Gmod.z,'b-',1:length(Gmod.x),Gmod.r,'k-')            subplot(2,2,3)            plot(1:td1,Gnorm.x,'ro',1:td1,Gnorm.y,'go',1:td1,Gnorm.z,'bo')            Gdir = Gnorm.x + i*Gnorm.y;            subplot(2,2,4)            plot(real(Gdir),imag(Gdir),'bo')            [X,Y] = fSchmidt(Gnorm.x,Gnorm.y,Gnorm.z);            latitude.theta = pi/180*[30:30:150 179];            latitude.phi = linspace(0,2*pi,100);            [latitude.phi,latitude.theta] = ndgrid(latitude.phi,latitude.theta);            latitude.z = cos(latitude.theta);            latitude.x = sin(latitude.theta).*cos(latitude.phi);            latitude.y = sin(latitude.theta).*sin(latitude.phi);            [latitude.X,latitude.Y] = fSchmidt(latitude.x,latitude.y,latitude.z);            longitude.theta = pi/180*linspace(30,180,100);            longitude.phi = pi/180*[30:30:360];            [longitude.theta,longitude.phi] = ndgrid(longitude.theta,longitude.phi);            longitude.z = cos(longitude.theta);            longitude.x = sin(longitude.theta).*cos(longitude.phi);            longitude.y = sin(longitude.theta).*sin(longitude.phi);            [longitude.X,longitude.Y] = fSchmidt(longitude.x,longitude.y,longitude.z);            subplot(2,2,4)                    plot(X,Y,'ko')            hold on            plot(latitude.X,latitude.Y,'b-')            plot(longitude.X,longitude.Y,'b-')            axis equal                           if exist('td1end') == 0                td1end = td1;            end%%            eval(['load ' num2str(expno(nexp)) '/Images'])            figure(1), clf            for np = 1:length(Pixel.x)                               nj = Pixel.y(np);                    ni = Pixel.x(np)                        index_fit = td1start:min([length(G.r) td1end]);                        %Yin = squeeze(sum(sum(abs(Itd1(:,:,index_fit)),1),2));                        Yin = abs(squeeze(Itd1(ni,nj,index_fit)));                        Iplot = abs(squeeze(Itd1(:,:,(NMRacqu2s.td/Ndir+1))));                        Iplot = Iplot/max(max(Iplot));                        if Iplot(ni,nj)>thresh                            Xin = b(index_fit);                            Gnormpoints.x = Gnorm.x(index_fit); Gnormpoints.y = Gnorm.y(index_fit); Gnormpoints.z = Gnorm.z(index_fit);                                                        Yin_iso = zeros(NMRacqu2s.td/Ndir,1);                            Yin_aniso = zeros(NMRacqu2s.td/Ndir,1);                            Ndir_iso = 0; Ndir_aniso = 0;                            for ndir = 1:Ndir                                index_fit = ndir - Ndir + (Ndir:Ndir:length(G.r));                                index_fit = index_fit(1:length(index_fit));                                Yin = squeeze(Itd1(ni,nj,index_fit));                                Xin = b(index_fit);                                Gnormpoints.x = Gnorm.x(index_fit); Gnormpoints.y = Gnorm.y(index_fit); Gnormpoints.z = Gnorm.z(index_fit);                                if isnan(Gnormpoints.x(1)) == 1                                    if ndir == 1                                        Yin(1) = 0;                                    end                                    Yin_iso = Yin_iso+Yin;                                    Xin_iso = Xin;                                    Ndir_iso = Ndir_iso+1;                                else                                    Yin_aniso = Yin_aniso+Yin;                                    Xin_aniso = Xin;                                    Ndir_aniso = Ndir_aniso+1;                                end                            end                            Yin_aniso = abs(Yin_aniso)/Ndir_aniso;                            Yin_iso = abs(Yin_iso)/Ndir_iso;                            Yin_iso(1) = Yin_iso(1)*Ndir_iso/(Ndir_iso-1);                            if Ndir_iso==1, Yin_iso(1) = Yin_aniso(1); end                                                        S0_gamma = Images.S0_gamma(ni,nj);                            mD_gamma = Images.mD_gamma(ni,nj);                            mu2_aniso = Images.mu2_aniso(ni,nj);                            mu2_iso = Images.mu2_iso(ni,nj);                            mu2_FA = 4/45*((Images.lambda.z(ni,nj)-Images.lambda.x(ni,nj))^2+...                                (Images.lambda.y(ni,nj)-Images.lambda.z(ni,nj))...                                *(Images.lambda.y(ni,nj)-Images.lambda.x(ni,nj)));                            FA = Images.FA(ni,nj);                            uFA = Images.uFA(ni,nj);                            K = Images.K(ni,nj);                                                        Ycalc_mD = 1*exp(-mD_gamma*Xin_iso);                            Ycalc_aniso = fexpgamma([1 mD_gamma sqrt(mu2_aniso)],Xin_aniso);                            Ycalc_iso = fexpgamma([1 mD_gamma sqrt(mu2_iso)],Xin_iso);                            Ycalc_FA = fexpgamma([1 mD_gamma sqrt(mu2_FA)],Xin_iso);                                                                subplot(2,2,np)                            semilogy(Xin_iso*mD_gamma,Ycalc_mD,'k--','LineWidth',2*lw)                            hold on                            %semilogy(Xin_iso*mD_gamma,Ycalc_FA,'k--','LineWidth',2*lw)                            semilogy(Xin_aniso*mD_gamma,Ycalc_aniso,'b-','LineWidth',lw)                            semilogy(Xin_aniso*mD_gamma,Yin_aniso/S0_gamma,...                                'bo','LineWidth',lw,'MarkerSize',8)                            semilogy(Xin_iso*mD_gamma,Yin_iso/S0_gamma,'ro',...                                'LineWidth',lw,'MarkerSize',5,'MarkerFaceColor','r')                            semilogy(Xin_iso*mD_gamma,Ycalc_iso,'r-','LineWidth',lw)                            hold off                            set(gca,'XLim',log(10)*[-.1 1.1],'YLim',[.09 1.1])                                                        set(gca,'Box','off','TickDir','out','TickLength',tl*[1 1],...                                'FontSize',fs,'LineWidth',lw)                            title(['FA=' num2str(FA,2) ' uFA=' num2str(uFA,2)])                                                end                end                   end    end    eval(['print ' num2str(expno(nexp)) '/Ivsb_pixel -loose -depsc'])    cd ..endcd(wd)