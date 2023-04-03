clear allwd = cd;%DataDir = '/opt/topspin2/data/DT/nmr';DataDir = '/Users/daniel/NMRdata/AVII500/DT';%DataDir = '/Users/daniel/Documents/Spaces/Presentations';%DataDir = '/Users/daniel/Dropbox';%ExpNam = {'DDeltaMap'}; expno = 65:68;%ExpNam = {'AOToct8'}; expno = 17:18; expno = 11;ExpNam = {'AOToct9'}; expno = 17; expno = 11;%ExpNam = {'C14E5_6'}; expno = 20; expno = 11; expno = 17;%ExpNam = {'AOToct8'; 'AOToct9'}; expno = 17;ExpNam = {'C14E5_7'}; expno = 11;fs = 15;cd(DataDir)%GetExpnamsDmax = 3e-9;for ndir = 1:length(ExpNam)    ExpNam{ndir}    cd(ExpNam{ndir})    if exist('expno') == 0        GetExpnos    end    for nexp = 1:length(expno)        eval(['load ' num2str(expno(nexp)) '/NMRacqus'])        eval(['load ' num2str(expno(nexp)) '/Images'])                    EllipseScale = .3;                figure(4), clf        axes('position',[1-8/11 0 8/11 1])        view(0,90)        hold on        nSQ = 0;        SQ.MD = zeros(numel(Images.lambda1),1);        SQ.DDelta = zeros(numel(Images.lambda1),1);        SQ.lambda1 = zeros(numel(Images.lambda1),1);        SQ.lambda2 = zeros(numel(Images.lambda1),1);        SQ.lambda3 = zeros(numel(Images.lambda1),1);        SQ.Ealpha = zeros(numel(Images.lambda1),1);        SQ.Ebeta = zeros(numel(Images.lambda1),1);        SQ.Egamma = zeros(numel(Images.lambda1),1);        SQ.x = zeros(numel(Images.lambda1),1);        SQ.y = zeros(numel(Images.lambda1),1);        SQ.z = zeros(numel(Images.lambda1),1);                dx = r.i(2) - r.i(1); dy = r.j(2) - r.j(1);        r.i = circshift(r.i,1);        r.j = circshift(r.j,0);        for ni = 1:numel(r.i)            for nj = 1:numel(r.j)                %ni = 3; nj = 3;                lambda.x = Images.lambda.x(ni,nj);                lambda.y = Images.lambda.y(ni,nj);                lambda.z = Images.lambda.z(ni,nj);                if all([lambda.x lambda.y lambda.z] > 0)                    A.alpha = Images.A.alpha(ni,nj);                    A.beta = Images.A.beta(ni,nj);                    A.gamma = Images.A.gamma(ni,nj);                    C.x = [1; 0; 0];                    C.y = [0; 1; 0];                    C.z = [0; 0; 1];                    R.gamma = [                    cos(A.gamma) sin(A.gamma) 0                    -sin(A.gamma) cos(A.gamma) 0                    0 0 1];                    R.beta = [                    cos(A.beta) 0 -sin(A.beta)                    0 1 0                    sin(A.beta) 0 cos(A.beta)];                    R.alpha = [                    cos(A.alpha) sin(A.alpha) 0                    -sin(A.alpha) cos(A.alpha) 0                    0 0 1];                    R.mat = R.gamma*R.beta*R.alpha;                    R.inv = inv(R.mat);                    gamma = .5;                    Nphi = 100;                    Ntheta = Nphi;                    phi = linspace(0,pi,Nphi);                    theta = linspace(0,2*pi,Ntheta);                    [phi,theta] = ndgrid(phi,theta);                    lambda1 = Images.lambda1(ni,nj);                    lambda2 = Images.lambda2(ni,nj);                    lambda3 = Images.lambda3(ni,nj);                    v1.x = Images.v1.x(ni,nj); v1.y = Images.v1.y(ni,nj); v1.z = Images.v1.z(ni,nj);                    v2.x = Images.v2.x(ni,nj); v2.y = Images.v2.y(ni,nj); v2.z = Images.v2.z(ni,nj);                    v3.x = Images.v3.x(ni,nj); v3.y = Images.v3.y(ni,nj); v3.z = Images.v3.z(ni,nj);                    v1.x = v3.y*v2.z - v3.z*v2.y;                    v1.y = v3.z*v2.x - v3.x*v2.z;                    v1.z = v3.x*v2.y - v3.y*v2.x;                                        cl = (lambda1-lambda2)/(lambda1+lambda2+lambda3);                    cp = 2*(lambda2-lambda3)/(lambda1+lambda2+lambda3);                    cs = 3*lambda3/(lambda1+lambda2+lambda3);                    if cl >= cp                        alpha = (1-cp)^gamma;                        beta = (1-cl)^gamma;                        q.x = sign(cos(phi)).*abs(cos(phi)).^beta;                        q.y = (-1)*sign(sin(theta)).*abs(sin(theta)).^alpha.*sign(sin(phi)).*abs(sin(phi)).^beta;                        q.z = sign(cos(theta)).*abs(cos(theta)).^alpha.*sign(sin(phi)).*abs(sin(phi)).^beta;                    elseif cl < cp                        alpha = (1-cl)^gamma;                        beta = (1-cp)^gamma;                        q.x = sign(cos(theta)).*abs(cos(theta)).^alpha.*sign(sin(phi)).*abs(sin(phi)).^beta;                        q.y = sign(sin(theta)).*abs(sin(theta)).^alpha.*sign(sin(phi)).*abs(sin(phi)).^beta;                        q.z = sign(cos(phi)).*abs(cos(phi)).^beta;                    end                    if lambda1 == lambda.x                                    if lambda2 == lambda.y                            DTsuperquad_PAS.x = -q.x;                            DTsuperquad_PAS.y = q.y;                            DTsuperquad_PAS.z = q.z;                        elseif lambda2 == lambda.z                            DTsuperquad_PAS.x = -q.x;                            DTsuperquad_PAS.z = -q.y;                            DTsuperquad_PAS.y = q.z;                        end                    elseif lambda1 == lambda.y                                    if lambda2 == lambda.x                            DTsuperquad_PAS.y = -q.x;                            DTsuperquad_PAS.x = q.y;                            DTsuperquad_PAS.z = -q.z;                        elseif lambda2 == lambda.z                            DTsuperquad_PAS.y = -q.x;                            DTsuperquad_PAS.z = q.y;                            DTsuperquad_PAS.x = q.z;                        end                    elseif lambda1 == lambda.z                                    if lambda2 == lambda.x                            DTsuperquad_PAS.z = -q.x;                            DTsuperquad_PAS.x = q.y;                            DTsuperquad_PAS.y = q.z;                        elseif lambda2 == lambda.y                            DTsuperquad_PAS.z = -q.x;                            DTsuperquad_PAS.y = q.y;                            DTsuperquad_PAS.x = -q.z;                        end                    end                    DTsuperquad_PAS.x = (lambda.x/Dmax)*DTsuperquad_PAS.x;                    DTsuperquad_PAS.y = (lambda.y/Dmax)*DTsuperquad_PAS.y;                    DTsuperquad_PAS.z = (lambda.z/Dmax)*DTsuperquad_PAS.z;                    DTsuperquad_LF.x = R.mat(1,1)*DTsuperquad_PAS.x + R.mat(1,2)*DTsuperquad_PAS.y + R.mat(1,3)*DTsuperquad_PAS.z;                    DTsuperquad_LF.y = R.mat(2,1)*DTsuperquad_PAS.x + R.mat(2,2)*DTsuperquad_PAS.y + R.mat(2,3)*DTsuperquad_PAS.z;                    DTsuperquad_LF.z = R.mat(3,1)*DTsuperquad_PAS.x + R.mat(3,2)*DTsuperquad_PAS.y + R.mat(3,3)*DTsuperquad_PAS.z;                    X = EllipseScale*DTsuperquad_LF.x + 1e3*r.i(ni);                    Y = EllipseScale*DTsuperquad_LF.y + 1e3*r.j(nj);                    Z = EllipseScale*DTsuperquad_LF.z;            %         S = [30,30];            %         h = surfl(X,Y,Z,S);            %         hold on                    h = surfl(X,Y,Z);%                     X = [0; EllipseScale*v1.x] + 1e3*r.i(ni);%                     Y = [0; EllipseScale*v1.y] + 1e3*r.j(nj);%                     Z = [0; EllipseScale*v1.z];%                     plot3(X,Y,Z,'b-')%                     X = [0; EllipseScale*v2.x] + 1e3*r.i(ni);%                     Y = [0; EllipseScale*v2.y] + 1e3*r.j(nj);%                     Z = [0; EllipseScale*v2.z];%                     plot3(X,Y,Z,'g-')%                     X = [0; EllipseScale*v3.x] + 1e3*r.i(ni);%                     Y = [0; EllipseScale*v3.y] + 1e3*r.j(nj);%                     Z = [0; EllipseScale*v3.z];%                     plot3(X,Y,Z,'r-')%                                         Ealpha = atan2(v1.y,v1.x);                    Ebeta = acos(v1.z);                     v1IF1.x = v1.x*cos(Ealpha) + v1.y*sin(Ealpha);                    v1IF1.y = v1.y*cos(Ealpha) - v1.x*sin(Ealpha);                    v1IF1.z = v1.z;                    v2IF1.x = v2.x*cos(Ealpha) + v2.y*sin(Ealpha);                    v2IF1.y = v2.y*cos(Ealpha) - v2.x*sin(Ealpha);                    v2IF1.z = v2.z;                    v3IF1.x = v3.x*cos(Ealpha) + v3.y*sin(Ealpha);                    v3IF1.y = v3.y*cos(Ealpha) - v3.x*sin(Ealpha);                    v3IF1.z = v3.z;                    v1IF2.x = v1IF1.x*cos(Ebeta) - v1IF1.z*sin(Ebeta);                    v1IF2.y = v1IF1.y;                    v1IF2.z = v1IF1.z*cos(Ebeta) + v1IF1.x*sin(Ebeta);                    v2IF2.x = v2IF1.x*cos(Ebeta) - v2IF1.z*sin(Ebeta);                    v2IF2.y = v2IF1.y;                    v2IF2.z = v2IF1.z*cos(Ebeta) + v2IF1.x*sin(Ebeta);                    v3IF2.x = v3IF1.x*cos(Ebeta) - v3IF1.z*sin(Ebeta);                    v3IF2.y = v3IF1.y;                    v3IF2.z = v3IF1.z*cos(Ebeta) + v3IF1.x*sin(Ebeta);                    Egamma = atan2(v3IF2.y,v3IF2.x);% %                     X = [0; EllipseScale*v1IF1.x] + 1e3*r.i(ni);%                     Y = [0; EllipseScale*v1IF1.y] + 1e3*r.j(nj);%                     Z = [0; EllipseScale*v1IF1.z];%                     plot3(X,Y,Z,'k-')%                     X = [0; EllipseScale*v1IF2.x] + 1e3*r.i(ni);%                     Y = [0; EllipseScale*v1IF2.y] + 1e3*r.j(nj);%                     Z = [0; EllipseScale*v1IF2.z];%                     plot3(X,Y,Z,'b--','LineWidth',2)%                     X = [0; 0] + 1e3*r.i(ni);%                     Y = [0; 0] + 1e3*r.j(nj);%                     Z = [0; 1.1*EllipseScale];%                     plot3(X,Y,Z,'k-','LineWidth',1)%                     X = [0; 1.1*EllipseScale] + 1e3*r.i(ni);%                     Y = [0; 0] + 1e3*r.j(nj);%                     Z = [0; 0];%                     plot3(X,Y,Z,'k-','LineWidth',1)%                     X = [0; 0] + 1e3*r.i(ni);%                     Y = [0; 1.1*EllipseScale] + 1e3*r.j(nj);%                     Z = [0; 0];%                     plot3(X,Y,Z,'k-','LineWidth',1)%                     X = [0; EllipseScale*v3IF2.x] + 1e3*r.i(ni);%                     Y = [0; EllipseScale*v3IF2.y] + 1e3*r.j(nj);%                     Z = [0; EllipseScale*v3IF2.z];%                     plot3(X,Y,Z,'r--','LineWidth',2)%                     X = [0; EllipseScale*v2IF2.x] + 1e3*r.i(ni);%                     Y = [0; EllipseScale*v2IF2.y] + 1e3*r.j(nj);%                     Z = [0; EllipseScale*v2IF2.z];%                     plot3(X,Y,Z,'g--','LineWidth',2)%                                         nSQ = nSQ+1;                    DDelta = Images.DDelta(ni,nj);                    MD = Images.MD_erf(ni,nj);                    SQ.MD(nSQ,1) = MD/1e-9;                    SQ.DDelta(nSQ,1) = DDelta;                    SQ.lambda1(nSQ,1) = lambda1/1e-9;                    SQ.lambda2(nSQ,1) = lambda2/1e-9;                    SQ.lambda3(nSQ,1) = lambda3/1e-9;                    SQ.Ealpha(nSQ,1) = Ealpha/pi*180;                    SQ.Ebeta(nSQ,1) = Ebeta/pi*180;                    SQ.Egamma(nSQ,1) = Egamma/pi*180;                    SQ.x(nSQ,1) = 1e3*r.i(ni);                    SQ.y(nSQ,1) = 1e3*r.j(nj);                    SQ.z(nSQ,1) = 0;                    %                     3*[lambda1 lambda2 lambda3]/(lambda1+lambda2+lambda3)%                     (3*[lambda1 lambda2 lambda3]/(lambda1+lambda2+lambda3) - 1)/2/DDelta                                    end                SQ.N = nSQ;                                %return            pause(.1)            end        end        axis tight        axis equal        %alpha(0.3)        shading flat        colormap(gray(256))        %colormap(copper(64))        %view(0,60)        xlabel('x'), ylabel('y')        set(gca,'XLim',1e3*[min(r.i)-dx/2 max(r.i)+dx/2])        set(gca,'YLim',1e3*[min(r.j)-dy/2 max(r.j)+dy/2])        %set(gca,'ZLim',1e3*[min(r.z)-dz/2 max(r.z)+dz/2])        axis off                fid = fopen([num2str(expno(nexp)) '/SQdat.txt'], 'w');        N = SQ.N;        fprintf(fid, '%8.0i,\n', N);        format = '%8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f,\n';        for n = 1:N            fprintf(fid,format,[SQ.MD(n) SQ.DDelta(n) SQ.lambda1(n) SQ.lambda2(n) SQ.lambda3(n) SQ.Ealpha(n) SQ.Ebeta(n) SQ.Egamma(n) SQ.x(n) SQ.y(n) SQ.z(n)]);        end        fclose(fid);                 end    cd ..endcd(wd)