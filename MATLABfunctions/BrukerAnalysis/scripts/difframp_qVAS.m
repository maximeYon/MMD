clear all

wd = cd;
DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
% ExpNam = 'AOT_1H2H'; expno = 87;
% 
% %ConeAngleUnit = pi/180*[54.7356*[.5 1] 90 180-54.7356*[1 .5 0]];
% ConeAngleUnit = pi/180*[54.7356*[1/3 2/3 1] 75 90 105 180-54.7356*[1 2/3 1/3 0]];
% %ConeAngleUnit = [ConeAngleUnit ConeAngleUnit+pi]
% %ConeAngleUnit = 2*pi*linspace(0,1,25); ConeAngleUnit = ConeAngleUnit(1:end-1);
% NConeCycles = 8;
% ConeAngleStep = 2*pi*((1:NConeCycles)-1);
% [ConeAngleUnit, ConeAngleStep] = ndgrid(ConeAngleUnit,ConeAngleStep);
% NG = numel(ConeAngleUnit);
% ConeAngle = reshape(ConeAngleUnit,NG,1) + reshape(ConeAngleStep,NG,1);
% 
% NG = 8;
% ConeAngle = pi/180*[0*54.7356]*ones(NG,1);
% %Grel = sqrt(linspace(.02,1,NG))';
% Grel = linspace(.02,1,NG)';
% %Grel = linspace(.999,1,NG)';
% 
% load UniformDistSphereRepulsionN7
% %theta = 0*pi/2; phi = 0*pi/2;
% Ndir = length(theta);
%  phi = pi/2*ones(Ndir,1);
%  theta = -pi/2 + pi/2*linspace(0,1,Ndir)';
% 
% [Gindx,Dindx] = ndgrid(1:NG,1:Ndir);
% td1 = numel(Gindx);
% 
% Gindx = reshape(Gindx,td1,1);
% Dindx = reshape(Dindx,td1,1);
% G.a = 1/sqrt(2)*sin(ConeAngle(Gindx)).*Grel(Gindx);
% G.b = 1/sqrt(2)*sin(ConeAngle(Gindx)).*Grel(Gindx);
% G.c = cos(ConeAngle(Gindx)).*Grel(Gindx);
% 
% alpha = phi(Dindx);
% beta = theta(Dindx);
% gamma = zeros(size(Dindx));

ExpNam = 'C14E5_5'; expno = 29; 
ExpNam = 'DDeltaMap'; expno = [3]; expno = 41; expno = 64;
ExpNam = 'AOToct8'; expno = 11; 
ExpNam = 'C14E5_temp1'; expno = 3;
ExpNam = 'Saupe_SliceTest'; expno = 4;
ExpNam = 'AOToct_Eq2'; expno = 6;
ExpNam = 'RARE2DAxSym_EddyTest'; expno = 29;
ExpNam = 'AOToct_Eq3'; expno = 6;
ExpNam = 'AOToct_Eq4'; expno = 6;

%relative b-values; 0 - 1

NbTtrace = 4;
brelmin = .01;
brelminlog = .1;
brelmaxiso = .4;
%brel = linspace(brelmin,1,NbTtrace)'; %linear spacing b
brel = logspace(log10(brelminlog),0,NbTtrace)'; brel(1) = brelmin; %geometric spacing b

%b-tensor anisotropy; -0.5 - 1

NbTDelta = 4; %NbTDelta = 4, 7, 10, 13,...
bTDelta = linspace(.9999,-.4999,NbTDelta)';
%bTDelta = [0]'; NbTDelta = length(bTDelta); %Sphere
%bTDelta = [1]'; NbTDelta = length(bTDelta); %Stick
%bTDelta = [0 1]'; NbTDelta = length(bTDelta); %Sphere and stick
%bTDelta = [1 0 -.5]'; NbTDelta = length(bTDelta); %Stick, sphere, and disk
%bTDelta = [1 1 1 .8]'; NbTDelta = length(bTDelta);


%directions

NbTdir = 17;
%electrostatic repulsion directions; NbTdir = 5, 6, ..., 97
%eval(['load UniformDistSphereRepulsionN' num2str(NbTdir)])
%alt algorithm; NbTdir = 10, 11, ..., 1000
eval(['load UDSRTriN' num2str(NbTdir)]) , phi = UDSR.phi; theta = UDSR.theta;
%single direction
%theta = 0*pi/2; phi = 0*pi/2; NbTdir = length(theta);


%make 3D grid of gradients, bTDelta, and directions

[trace_indx,bTDelta_indx,dir_indx] = ndgrid(1:NbTtrace,1:NbTDelta,1:NbTdir);
td1 = numel(trace_indx);

trace_indx = reshape(trace_indx,td1,1);
bTDelta_indx = reshape(bTDelta_indx,td1,1);
dir_indx = reshape(dir_indx,td1,1);

alpha = phi(dir_indx);
beta = theta(dir_indx);
gamma = zeros(size(dir_indx));
%figure(1), clf, plot(alpha,beta,'o'), return

Grel = sqrt(brel);
ConeAngle = acos(sqrt((bTDelta*2+1)/3));


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
 
% G.a = 1/sqrt(2)*sin(ConeAngle(bTDelta_indx)).*Grel(trace_indx);
% G.b = 1/sqrt(2)*sin(ConeAngle(bTDelta_indx)).*Grel(trace_indx);
% G.a = sin(ConeAngle(bTDelta_indx)).*Grel(trace_indx);
% G.b = sin(ConeAngle(bTDelta_indx)).*Grel(trace_indx);
% G.c = cos(ConeAngle(bTDelta_indx)).*Grel(trace_indx);
% G.a = sin(ConeAngle(bTDelta_indx)).*Grel(trace_indx)...
%     .*(sqrt(brelmaxiso) + (1-sqrt(brelmaxiso)).*abs(bTDelta(bTDelta_indx)));
% G.b = sin(ConeAngle(bTDelta_indx)).*Grel(trace_indx)...
%     .*(sqrt(brelmaxiso) + (1-sqrt(brelmaxiso)).*abs(bTDelta(bTDelta_indx)));
% G.c = cos(ConeAngle(bTDelta_indx)).*Grel(trace_indx)...
%     .*(sqrt(brelmaxiso) + (1-sqrt(brelmaxiso)).*abs(bTDelta(bTDelta_indx)));
G.a = sin(ConeAngle(bTDelta_indx)).*Grel(trace_indx)...
    .*(sqrt(brelmaxiso) + (1-sqrt(brelmaxiso)).*abs(bTDelta(bTDelta_indx)).^2);
G.b = sin(ConeAngle(bTDelta_indx)).*Grel(trace_indx)...
    .*(sqrt(brelmaxiso) + (1-sqrt(brelmaxiso)).*abs(bTDelta(bTDelta_indx)).^2);
G.c = cos(ConeAngle(bTDelta_indx)).*Grel(trace_indx)...
    .*(sqrt(brelmaxiso) + (1-sqrt(brelmaxiso)).*abs(bTDelta(bTDelta_indx)).^2);

rax = R.xx.*G.a;
ray = R.xy.*G.a;
raz = R.xz.*G.a;
rbx = R.yx.*G.b;
rby = R.yy.*G.b;
rbz = R.yz.*G.b;
rcx = R.zx.*G.c;
rcy = R.zy.*G.c;
rcz = R.zz.*G.c;

rtot = 0*sqrt((rax+rbx+rcx).^2+(ray+rby+rcy).^2+(raz+rbz+rcz).^2);

figure(1), clf
subplot(1,2,1)
plot(1:td1,rax,'ro',1:td1,rbx,'rs',1:td1,rcx,'rx',...
    1:td1,ray,'go',1:td1,rby,'gs',1:td1,rcy,'gx',...
    1:td1,raz,'bo',1:td1,rbz,'bs',1:td1,rcz,'bx',...
    1:td1,rtot,'k-')

ylabel('g ramps')

title(['td1=' num2str(td1)])

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

subplot(1,2,2)        
plot(X,Y,'ko')
hold on
plot(latitude.X,latitude.Y,'b-')
plot(longitude.X,longitude.Y,'b-')
axis tight equal off
%return

RampMat = struct('ax',rax,'ay',ray,'az',raz,'bx',rbx,'by',rby,'bz',rbz,'cx',rcx,'cy',rcy,'cz',rcz);
dir = struct('x',Gnorm.x,'y',Gnorm.y,'z',Gnorm.z,'theta',theta,'phi',phi);
DiffRamp = struct('GridSamp',1,'NbTtrace',NbTtrace,'NbTDelta',NbTDelta,'NbTdir',NbTdir,...
'mat',RampMat,'dir',dir);

fpath = [DataDir '/' ExpNam '/' num2str(expno)];
eval(['save ' fpath '/DiffRamp DiffRamp'])

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
{'##TITLE= /opt/topspin2/data/DT/nmr/C10E3_DTI/6/diff_ramp'};
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
