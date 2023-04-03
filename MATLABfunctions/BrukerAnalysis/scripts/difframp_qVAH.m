clear all

wd = cd;
DataDir = '/opt/topspin2/data/DT/nmr';
DataDir = '/Users/daniel/NMRdata/AVII500/DT';

%ExpNam = 'Yeast_tripgse'; expno = 64;
ExpNam = 'dummy'; expno = 1;

NDeltab = 4;
Deltab = linspace(.96,-.48,NDeltab)';
%Deltab = linspace(1,-.5,NDeltab)';
%NDeltab = 2; Deltab = [0 .96]';
%Deltab = [0.96 0 -.48]'; NDeltab = 3;
ConeAngle = acos(sqrt((Deltab*2+1)/3));

NG = 8;
Grel = sqrt(linspace(.2^2,1,NG))';
%Grel = sqrt(logspace(log10(0.02),0,NG))';
%Grel = sqrt(logspace(log10(0.01),0,NG))';
%Grel = sqrt(logspace(log10(0.005),0,NG))';
%Grel = sqrt(logspace(log10(0.001),0,NG))';
minG = 0;

Ndir = 17;
eval(['load UniformDistSphereRepulsionN' num2str(Ndir)])
%theta = 0*pi/2; phi = 0*pi/2; Ndir = length(theta);
%phi = 0*pi/2*ones(Ndir,1);
%theta = -0*pi/2 + pi/2*linspace(0,1,Ndir)';

[G_indx,Deltab_indx,Dir_indx] = ndgrid(1:NG,1:NDeltab,1:Ndir);
td1 = numel(G_indx);

G_indx = reshape(G_indx,td1,1);
Deltab_indx = reshape(Deltab_indx,td1,1);
Dir_indx = reshape(Dir_indx,td1,1);

gamma = phi(Dir_indx);
beta = theta(Dir_indx);
alpha = zeros(size(Dir_indx));
%figure(1), clf, plot(alpha,beta,'o'), return

for ntd1 = 1:td1
    A.gamma = gamma(ntd1);
    A.beta = beta(ntd1);
    A.alpha = alpha(ntd1);
    
    R.gamma = [
        cos(A.gamma) sin(A.gamma) 0
        -sin(A.gamma) cos(A.gamma) 0
        0 0 1];

    R.beta = [
        cos(A.beta) 0 -sin(A.beta)
        0 1 0
        sin(A.beta) 0 cos(A.beta)];

    R.alpha = [
        cos(A.alpha) sin(A.alpha) 0
        -sin(A.alpha) cos(A.alpha) 0
        0 0 1];

    R.mat = R.gamma*R.beta*R.alpha;
    
    R.xx(ntd1,1) = R.mat(1,1);
    R.xy(ntd1,1) = R.mat(1,2);
    R.xz(ntd1,1) = R.mat(1,3);
    R.yx(ntd1,1) = R.mat(2,1);
    R.yy(ntd1,1) = R.mat(2,2);
    R.yz(ntd1,1) = R.mat(2,3);
    R.zx(ntd1,1) = R.mat(3,1);
    R.zy(ntd1,1) = R.mat(3,2);
    R.zz(ntd1,1) = R.mat(3,3);
end
                
G.ax = cos(0).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
G.ay = sin(0).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
G.az = cos(ConeAngle(Deltab_indx)).*Grel(G_indx);
G.bx = cos(2*pi/3).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
G.by = sin(2*pi/3).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
G.bz = cos(ConeAngle(Deltab_indx)).*Grel(G_indx);
G.cx = cos(4*pi/3).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
G.cy = sin(4*pi/3).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
G.cz = cos(ConeAngle(Deltab_indx)).*Grel(G_indx);

rax = R.xx.*G.ax + R.xy.*G.ay + R.xz.*G.az;
ray = R.yx.*G.ax + R.yy.*G.ay + R.yz.*G.az;
raz = R.zx.*G.ax + R.zy.*G.ay + R.zz.*G.az;
rbx = R.xx.*G.bx + R.xy.*G.by + R.xz.*G.bz;
rby = R.yx.*G.bx + R.yy.*G.by + R.yz.*G.bz;
rbz = R.zx.*G.bx + R.zy.*G.by + R.zz.*G.bz;
rcx = R.xx.*G.cx + R.xy.*G.cy + R.xz.*G.cz;
rcy = R.yx.*G.cx + R.yy.*G.cy + R.yz.*G.cz;
rcz = R.zx.*G.cx + R.zy.*G.cy + R.zz.*G.cz;

rax(abs(rax)<minG) = minG;
ray(abs(ray)<minG) = minG;
raz(abs(raz)<minG) = minG;
rbx(abs(rbx)<minG) = minG;
rby(abs(rby)<minG) = minG;
rbz(abs(rbz)<minG) = minG;
rcx(abs(rcx)<minG) = minG;
rcy(abs(rcy)<minG) = minG;
rcz(abs(rcz)<minG) = minG;

symv.x = rax+rbx+rcx;
symv.y = ray+rby+rcy;
symv.z = raz+rbz+rcz;


rtot = sqrt((rax+rbx+rcx).^2+(ray+rby+rcy).^2+(raz+rbz+rcz).^2);

figure(1), clf
subplot(2,2,1)
plot(1:td1,rax,'ro',1:td1,rbx,'rs',1:td1,rcx,'rx',...
    1:td1,ray,'go',1:td1,rby,'gs',1:td1,rcy,'gx',...
    1:td1,raz,'bo',1:td1,rbz,'bs',1:td1,rcz,'bx')
ylabel('g ramps')

title(['td1=' num2str(td1)])

bmat.xx = rax.*rax + rbx.*rbx + rcx.*rcx;
bmat.yy = ray.*ray + rby.*rby + rcy.*rcy;
bmat.zz = raz.*raz + rbz.*rbz + rcz.*rcz;
bmat.xy = rax.*ray + rbx.*rby + rcx.*rcy;
bmat.xz = rax.*raz + rbx.*rbz + rcx.*rcz;
bmat.yz = ray.*raz + rby.*rbz + rcy.*rcz;

btot = bmat.xx + bmat.yy + bmat.zz;

lambda1 = zeros(1,td1);
lambda2 = zeros(1,td1);
lambda3 = zeros(1,td1);
for ntd1 = 1:td1
    lambdas = eig([bmat.xx(ntd1) bmat.xy(ntd1) bmat.xz(ntd1)
    bmat.xy(ntd1) bmat.yy(ntd1) bmat.yz(ntd1)
    bmat.xz(ntd1) bmat.yz(ntd1) bmat.zz(ntd1)]);
    lambda.tot(1,ntd1) = sum(lambdas);
    Dlambdas = abs(lambdas-lambda.tot(1,ntd1)/3);
    [dummy,indx] = sort(Dlambdas,'descend');
    lambda.c(1,ntd1) = lambdas(indx(1));
    lambda.b(1,ntd1) = lambdas(indx(2));
    lambda.a(1,ntd1) = lambdas(indx(3));                        
end

lambda.Delta = (lambda.c - (lambda.b+lambda.a)/2)./lambda.tot;
lambda.eta = (lambda.b - lambda.a)./lambda.tot;
%%
subplot(4,2,5)
plot(1:td1,lambda.tot)
ylabel('b')
subplot(4,2,7)
plot(1:td1,lambda.Delta)
ylabel('\Delta_b')
xlabel('td1')
%%
Gnorm.x = sin(theta).*cos(phi);
Gnorm.y = sin(theta).*sin(phi);
Gnorm.z = cos(theta);

[X,Y] = fSchmidt(Gnorm.x,Gnorm.y,Gnorm.z);
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

subplot(2,2,2)        
plot(X,Y,'ko')
hold on
plot(latitude.X,latitude.Y,'b-')
plot(longitude.X,longitude.Y,'b-')
axis tight equal off
title(['Ndir=' num2str(Ndir)])

subplot(2,2,4)
indx = NG+(0:(Ndir-1))*NG*NDeltab;
plot3([zeros(1,numel(indx)); symv.x(indx)'],[zeros(1,numel(indx)); symv.y(indx)'],[zeros(1,numel(indx)); symv.z(indx)'],'-o')
axis tight square
xlabel('x')
ylabel('y')
zlabel('z')
%return

%return

DRdat = struct('GridSamp',1,'Ndir',Ndir,'NDeltab',NDeltab);

fpath = [DataDir '/' ExpNam '/' num2str(expno)];
eval(['save ' fpath '/DRdat DRdat'])

fname = 'rax';
fidax = fopen([fpath '/' fname],'w');
fname = 'rax.txt';
fidaxtxt = fopen([fpath '/' fname],'w');
fname = 'ray';
fiday = fopen([fpath '/' fname],'w');
fname = 'ray.txt';
fidaytxt = fopen([fpath '/' fname],'w');
fname = 'raz';
fidaz = fopen([fpath '/' fname],'w');
fname = 'raz.txt';
fidaztxt = fopen([fpath '/' fname],'w');
fname = 'rbx';
fidbx = fopen([fpath '/' fname],'w');
fname = 'rbx.txt';
fidbxtxt = fopen([fpath '/' fname],'w');
fname = 'rby';
fidby = fopen([fpath '/' fname],'w');
fname = 'rby.txt';
fidbytxt = fopen([fpath '/' fname],'w');
fname = 'rbz';
fidbz = fopen([fpath '/' fname],'w');
fname = 'rbz.txt';
fidbztxt = fopen([fpath '/' fname],'w');
fname = 'rcx';
fidcx = fopen([fpath '/' fname],'w');
fname = 'rcx.txt';
fidcxtxt = fopen([fpath '/' fname],'w');
fname = 'rcy';
fidcy = fopen([fpath '/' fname],'w');
fname = 'rcy.txt';
fidcytxt = fopen([fpath '/' fname],'w');
fname = 'rcz';
fidcz = fopen([fpath '/' fname],'w');
fname = 'rcz.txt';
fidcztxt = fopen([fpath '/' fname],'w');

text.header = {
{['##TITLE= ' fpath '/diff_ramp']};
{'##JCAMP-DX= 5.00 Bruker JCAMP library'}
{'##DATA TYPE= Shape Data'}
{'##ORIGIN= Bruker Analytik GmbH'}
{'##OWNER= <nmrsu>'}
{'##DATE= 07/02/27'}
{'##TIME= 08:11:48'}
{'##MINX= 0'}
{'##MAXX= 1'}
{'##MINY= 0'}
{'##MAXY= 1'}
{'##$SHAPE_EXMODE= Gradient'}
{'##$SHAPE_TOTROT= 0'}
{'##$SHAPE_BWFAC= 0'}
{'##$SHAPE_INTEGFAC= 0'}
{'##$SHAPE_MODE= 0'}};

[Nlines, dummy] = size(text.header);

for nline = 1:Nlines
    fprintf(fidax,'%s\n',text.header{nline}{1});
    fprintf(fiday,'%s\n',text.header{nline}{1});
    fprintf(fidaz,'%s\n',text.header{nline}{1});
    fprintf(fidbx,'%s\n',text.header{nline}{1});
    fprintf(fidby,'%s\n',text.header{nline}{1});
    fprintf(fidbz,'%s\n',text.header{nline}{1});
    fprintf(fidcx,'%s\n',text.header{nline}{1});
    fprintf(fidcy,'%s\n',text.header{nline}{1});
    fprintf(fidcz,'%s\n',text.header{nline}{1});
end


fprintf(fidax,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidax,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fiday,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fiday,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidaz,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidaz,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidbx,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidbx,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidby,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidby,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidbz,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidbz,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidcx,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidcx,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidcy,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidcy,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidcz,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidcz,'%s\n',['##XYDATA= (X++(Y..Y))']);

for nG = 1:length(rax)
    fprintf(fidax,'%f\n',rax(nG));
    fprintf(fidaxtxt,'%f\n',rax(nG));
    fprintf(fiday,'%f\n',ray(nG));
    fprintf(fidaytxt,'%f\n',ray(nG));
    fprintf(fidaz,'%f\n',raz(nG));
    fprintf(fidaztxt,'%f\n',raz(nG));
    fprintf(fidbx,'%f\n',rbx(nG));
    fprintf(fidbxtxt,'%f\n',rbx(nG));
    fprintf(fidby,'%f\n',rby(nG));
    fprintf(fidbytxt,'%f\n',rby(nG));
    fprintf(fidbz,'%f\n',rbz(nG));
    fprintf(fidbztxt,'%f\n',rbz(nG));
    fprintf(fidcx,'%f\n',rcx(nG));
    fprintf(fidcxtxt,'%f\n',rcx(nG));
    fprintf(fidcy,'%f\n',rcy(nG));
    fprintf(fidcytxt,'%f\n',rcy(nG));
    fprintf(fidcz,'%f\n',rcz(nG));
    fprintf(fidcztxt,'%f\n',rcz(nG));
end


fprintf(fidax,'%s\n',['##END']);
fprintf(fiday,'%s\n',['##END']);
fprintf(fidaz,'%s\n',['##END']);
fprintf(fidbx,'%s\n',['##END']);
fprintf(fidby,'%s\n',['##END']);
fprintf(fidbz,'%s\n',['##END']);
fprintf(fidcx,'%s\n',['##END']);
fprintf(fidcy,'%s\n',['##END']);
fprintf(fidcz,'%s\n',['##END']);

fclose(fidax);
fclose(fiday);
fclose(fidaz);
fclose(fidaxtxt);
fclose(fidaytxt);
fclose(fidaztxt);
fclose(fidbx);
fclose(fidby);
fclose(fidbz);
fclose(fidbxtxt);
fclose(fidbytxt);
fclose(fidbztxt);
fclose(fidcx);
fclose(fidcy);
fclose(fidcz);
fclose(fidcxtxt);
fclose(fidcytxt);
fclose(fidcztxt);

td1