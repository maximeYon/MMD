clear all

wd = cd;

DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';

ExpNam = 'dummy'; expno = 1;
ExpNam = 'RandSampDR1R2'; expno = 140;
ExpNam = 'AOToct_temp10'; expno = 6;
ExpNam = 'AOToct8'; expno = 6;

%relative b-values; 0 - 1

Nb = 8;
brelmin = .002;
%brel = linspace(brelmin,1,Nb)'; %linear spacing b
brel = logspace(log10(brelmin),0,Nb)'; %geometric spacing b

%b-tensor anisotropy; -0.5 - 1

NDeltab = 4; %NDeltab = 4, 7, 10, 13,...
Deltab = linspace(1,-.5,NDeltab)';
%Deltab = [0]'; NDeltab = length(Deltab); %Sphere
%Deltab = [0 .96]'; NDeltab = length(Deltab); %Sphere and stick
%Deltab = [1 0 -.5]'; NDeltab = length(Deltab); %Stick, sphere, and disk

%directions

Ndir = 17;
%electrostatic repulsion directions; Ndir = 5, 6, ..., 97
%eval(['load UniformDistSphereRepulsionN' num2str(Ndir)])
%alt algorithm; Ndir = 10, 11, ..., 1000
eval(['load UDSRTriN' num2str(Ndir)]) , phi = UDSR.phi; theta = UDSR.theta;
%single direction
%theta = 1*pi/2; phi = 1*pi/2; Ndir = length(theta);


%make 3D grid of gradients, Deltab, and directions

[G_indx,Deltab_indx,Dir_indx] = ndgrid(1:Nb,1:NDeltab,1:Ndir);
td1 = numel(G_indx);

G_indx = reshape(G_indx,td1,1);
Deltab_indx = reshape(Deltab_indx,td1,1);
Dir_indx = reshape(Dir_indx,td1,1);

gamma = phi(Dir_indx);
beta = theta(Dir_indx);
alpha = zeros(size(Dir_indx));
%figure(1), clf, plot(alpha,beta,'o'), return

brel = brel(G_indx);
Deltab = Deltab(Deltab_indx);

bT.trace = brel;
bmin = brelmin;
bT.Delta = Deltab;
bT.par = (bT.trace-bmin)/3.*(1 + 2*bT.Delta) + bmin/3;
bT.perp = (bT.trace-bmin)/3.*(1 - bT.Delta) + bmin/3;
bT.trace = bT.par + 2*bT.perp;
bT.Delta = (3*bT.par./bT.trace - 1)/2;
bT.trace = bT.trace/max(bT.trace);
bT.par = bT.trace/3.*(1 + 2*bT.Delta);
bT.perp = bT.trace/3.*(1 - bT.Delta);
bT.s = 3*min([bT.par bT.perp],[],2);
bT.l = bT.par - bT.s/3;
bT.l(bT.l < 0) = 0;
bT.p = bT.trace - bT.s - bT.l;
brel = bT.trace;
Deltab = bT.Delta;

Grel = sqrt(brel);
ConeAngle = acos(sqrt((Deltab*2+1)/3));

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
                
% G.ax = cos(0).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
% G.ay = sin(0).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
% G.az = cos(ConeAngle(Deltab_indx)).*Grel(G_indx);
% G.bx = cos(2*pi/3).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
% G.by = sin(2*pi/3).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
% G.bz = cos(ConeAngle(Deltab_indx)).*Grel(G_indx);
% G.cx = cos(4*pi/3).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
% G.cy = sin(4*pi/3).*sin(ConeAngle(Deltab_indx)).*Grel(G_indx);
% G.cz = cos(ConeAngle(Deltab_indx)).*Grel(G_indx);

G.ax = cos(0).*sin(ConeAngle).*Grel;
G.ay = sin(0).*sin(ConeAngle).*Grel;
G.az = cos(ConeAngle).*Grel;
G.bx = cos(2*pi/3).*sin(ConeAngle).*Grel;
G.by = sin(2*pi/3).*sin(ConeAngle).*Grel;
G.bz = cos(ConeAngle).*Grel;
G.cx = cos(4*pi/3).*sin(ConeAngle).*Grel;
G.cy = sin(4*pi/3).*sin(ConeAngle).*Grel;
G.cz = cos(ConeAngle).*Grel;

rax = R.xx.*G.ax + R.xy.*G.ay + R.xz.*G.az;
ray = R.yx.*G.ax + R.yy.*G.ay + R.yz.*G.az;
raz = R.zx.*G.ax + R.zy.*G.ay + R.zz.*G.az;
rbx = R.xx.*G.bx + R.xy.*G.by + R.xz.*G.bz;
rby = R.yx.*G.bx + R.yy.*G.by + R.yz.*G.bz;
rbz = R.zx.*G.bx + R.zy.*G.by + R.zz.*G.bz;
rcx = R.xx.*G.cx + R.xy.*G.cy + R.xz.*G.cz;
rcy = R.yx.*G.cx + R.yy.*G.cy + R.yz.*G.cz;
rcz = R.zx.*G.cx + R.zy.*G.cy + R.zz.*G.cz;

symv.x = rax+rbx+rcx;
symv.y = ray+rby+rcy;
symv.z = raz+rbz+rcz;
% symv.norm = sqrt(symv.x.^2 + symv.y.^2 + symv.z.^2);
% symv.x = symv.x./symv.norm;
% symv.y = symv.y./symv.norm;
% symv.z = symv.z./symv.norm;
% symv.theta = acos(symv.z);
% symv.phi = atan2(symv.y,symv.x);

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
indx = Nb+(0:(Ndir-1))*Nb*NDeltab;
plot3([zeros(1,numel(indx)); symv.x(indx)'],[zeros(1,numel(indx)); symv.y(indx)'],[zeros(1,numel(indx)); symv.z(indx)'],'-o')
axis tight square
xlabel('x')
ylabel('y')
zlabel('z')

%save gradient ramps

DRdat = struct('GridSamp',1,'Nb',Nb,'NDeltab',NDeltab,'Ndir',Ndir);

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

cd(wd)