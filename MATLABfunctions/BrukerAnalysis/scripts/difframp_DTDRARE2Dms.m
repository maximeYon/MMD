clear all

wd = cd;

DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';

ExpNam = 'DTDRARE2Dms_test'; expno = 3;

%relative b-values; 0 - 1

sv.Nslices = 15;
fqmin = -70e3; fqmax = 70e3; 
fq = linspace(fqmin,fqmax,sv.Nslices)';

sv.NRDenc = 500;
%indxordered = [2:10:sv.NRDenc];
indxordered = [2:1:sv.NRDenc];

bmin = 0.01; bmax = 1;
vdT2min = 1e-3; vdT2max = 1e-3;
vdT1min = .01; vdT1max = .01;

%Acquisition protocol
sv.N = sv.Nslices*sv.NRDenc;
sv.bT.trace = bmin*(bmax/bmin).^rand(sv.N,1);
%bT.trace = bmin + (bmax-bmin)*rand(sv.N,1);
%sv.bT.trace = bmin + (bmax-bmin)*linspace(0,1,sv.NRDenc);
%sv.bT.trace = logspace(log10(bmin),log10(bmax),sv.NRDenc);
%sv.bT.trace = reshape(repmat(sv.bT.trace,[sv.Nslices 1]),[sv.N 1]);
sv.bT.Delta = -.5 + 1.5*rand(sv.N,1); %sv.bT.Delta = ones(sv.N,1);
sv.bT.theta = acos(2*rand(sv.N,1)-1);
sv.bT.phi = 2*pi*rand(sv.N,1);
sv.bT.gamma = 2*pi*rand(sv.N,1);
sv.vdT2 = vdT2min*(vdT2max/vdT2min).^rand(sv.N,1);
%sv.vdT2 = logspace(log10(vdT2min),log10(vdT2max),sv.NRDenc);
%sv.vdT2 = reshape(repmat(sv.vdT2,[sv.Nslices 1]),[sv.N 1]);
sv.vdT1 = vdT1min + (vdT1max-vdT1min)*rand(sv.N,1);

indx = [ceil(rand(1,ceil(.2*sv.N))*sv.N)];
sv.bT.Delta(indx) = 1;
sv.bT.trace(indx) = bmax;
sv.vdT2(indx) = .1*sv.vdT2(indx);
indx = [ceil(rand(1,ceil(.2*sv.N))*sv.N)];
sv.bT.Delta(indx) = 0;
sv.bT.trace(indx) = .3*sv.bT.trace(indx);
indx = [ceil(rand(1,ceil(.1*sv.N))*sv.N)];
sv.bT.Delta(indx) = -.5;
sv.bT.trace(indx) = bmax;
sv.vdT2(indx) = .1*sv.vdT2(indx);
indx = [1:(4*sv.Nslices) ceil(rand(1,ceil(.1*sv.N))*sv.N)];
sv.bT.Delta(indx) = 0;
sv.bT.trace(indx) = bmin;
sv.vdT2(indx) = vdT2min;
sv.vdT2(sv.vdT2<vdT2min) = vdT2min;

svunit.bT.trace = bmax*ones(sv.NRDenc,1);
svunit.bT.trace([1 5 9]) = bmin;
%svunit.bT.Delta = 1*ones(sv.NRDenc,1);
%svunit.bT.theta = 0*ones(sv.NRDenc,1);
%svunit.bT.phi = 0*ones(sv.NRDenc,1);
svunit.bT.Delta = [0 0 0 0 1 1 1 1 -.5 -.5 -.5 -.5];
svunit.bT.theta = pi/2*[0 0 1 1 0 0 1 1 0 0 1 1];
svunit.bT.phi = pi/2*[0 0 0 1 0 0 0 1 0 0 0 1];
svunit.bT.gamma = 0*ones(sv.NRDenc,1);
svunit.vdT1 = vdT1min*ones(sv.NRDenc,1);
svunit.vdT2 = vdT2min*ones(sv.NRDenc,1);

sv.fq = zeros(sv.N,1);
sv.sliceindx = zeros(sv.N,1);
for nRD = 1:sv.NRDenc
    indx = (nRD-1)*sv.Nslices + (1:sv.Nslices);
    indxrand = randperm(sv.Nslices);
    sv.sliceindx(indx) = indxrand;
    if any([nRD == [1]])
        sv.sliceindx(indx) = (1:sv.Nslices)';
        %sv.bT.trace(indx) = bmin;
        %sv.bT.Delta(indx) = 0;
        %sv.vdT2(indx) = vdT2min;
        %sv.vdT1(indx) = vdT1max;
    elseif any([nRD == [indxordered]])
        sv.sliceindx(indx) = (1:sv.Nslices)';
        %sv.bT.trace(indx) = bmin;
        %sv.bT.Delta(indx) = 0;
        %sv.vdT2(indx) = vdT2min;
        %sv.vdT1(indx) = vdT1max;
    elseif any([nRD == (indxordered-1)])
        sv.sliceindx(indx) = (1:sv.Nslices)';
        %sv.vdT1(indx) = vdT1max;
    end
%     sv.sliceindx(indx) = (1:sv.Nslices)';
%     sv.bT.trace(indx) = svunit.bT.trace(nRD);
%     sv.bT.Delta(indx) = svunit.bT.Delta(nRD);
%     sv.bT.theta(indx) = svunit.bT.theta(nRD);
%     sv.bT.phi(indx) = svunit.bT.phi(nRD);
%     sv.bT.gamma(indx) = svunit.bT.gamma(nRD);
%     sv.vdT2(indx) = svunit.vdT2(nRD);
%     sv.vdT1(indx) = svunit.vdT1(nRD);
end
sv.fq = fq(sv.sliceindx);

indx = sv.vdT2 > 20e-3;
sv.bT.trace(indx) = sv.bT.trace(indx)/10;
indx = sv.bT.trace < bmin;
sv.bT.trace(indx) = bmin;
%sv.bT.Delta(:) = 1;

% sv.bT.par = (sv.bT.trace-bmin)/3.*(1 + 2*sv.bT.Delta) + bmin/3;
% sv.bT.perp = (sv.bT.trace-bmin)/3.*(1 - sv.bT.Delta) + bmin/3;
% sv.bT.trace = sv.bT.par + 2*sv.bT.perp;
% sv.bT.Delta = (3*sv.bT.par./sv.bT.trace - 1)/2;
% sv.bT.trace = sv.bT.trace/max(sv.bT.trace);
% sv.bT.par = sv.bT.trace/3.*(1 + 2*sv.bT.Delta);
% sv.bT.perp = sv.bT.trace/3.*(1 - sv.bT.Delta);
% sv.bT.s = 3*min([sv.bT.par sv.bT.perp],[],2);
% sv.bT.l = sv.bT.par - sv.bT.s/3;
% sv.bT.l(sv.bT.l < 0) = 0;
% sv.bT.p = sv.bT.trace - sv.bT.s - sv.bT.l;

figure(2), clf
subplot(4,2,1)
plot(1:sv.N,sv.fq/1e3,'-')
ylabel('fq')
subplot(4,2,3)
semilogy(1:sv.N,sv.bT.trace,'.')
ylabel('b')
subplot(2,4,3)
semilogx(sv.bT.trace,sv.bT.Delta,'.')
xlabel('b'), ylabel('b_\Delta')
subplot(2,4,4)
loglog(sv.bT.trace,sv.vdT2,'.')
xlabel('b'), ylabel('vdT2')
subplot(2,2,3)
plot(sv.bT.theta,sv.bT.phi,'o')
xlabel('\theta'), ylabel('\phi')
subplot(4,2,6)
semilogy(1:sv.N,sv.vdT2,'.')
ylabel('vdT2')
subplot(4,2,8)
plot(1:sv.N,sv.vdT1,'.')
ylabel('vdT1')


sv.bT.dir.x = sin(sv.bT.theta).*cos(sv.bT.phi);
sv.bT.dir.y = sin(sv.bT.theta).*sin(sv.bT.phi);
sv.bT.dir.z = cos(sv.bT.theta);

Grel = sqrt(sv.bT.trace);
ConeAngle = acos(sqrt((sv.bT.Delta*2+1)/3));

td1 = sv.N;
alpha = sv.bT.phi;
beta = sv.bT.theta;
gamma = sv.bT.gamma;

%figure(1), clf, plot(alpha,beta,'o'), return

for ntd1 = 1:td1
    A.gamma = gamma(ntd1);
    A.beta = beta(ntd1);
    A.alpha = alpha(ntd1);
    
    R.gamma = [
        cos(A.gamma) -sin(A.gamma) 0
        sin(A.gamma) cos(A.gamma) 0
        0 0 1];

    R.beta = [
        cos(A.beta) 0 sin(A.beta)
        0 1 0
        -sin(A.beta) 0 cos(A.beta)];

    R.alpha = [
        cos(A.alpha) -sin(A.alpha) 0
        sin(A.alpha) cos(A.alpha) 0
        0 0 1];

    R.mat = R.alpha*R.beta*R.gamma;
    
    R.xa(ntd1,1) = R.mat(1,1);
    R.xb(ntd1,1) = R.mat(1,2);
    R.xc(ntd1,1) = R.mat(1,3);
    R.ya(ntd1,1) = R.mat(2,1);
    R.yb(ntd1,1) = R.mat(2,2);
    R.yc(ntd1,1) = R.mat(2,3);
    R.za(ntd1,1) = R.mat(3,1);
    R.zb(ntd1,1) = R.mat(3,2);
    R.zc(ntd1,1) = R.mat(3,3);
end
                
G.a = sin(ConeAngle).*Grel;
G.b = sin(ConeAngle).*Grel;
G.c = cos(ConeAngle).*Grel;

rxa = R.xa.*G.a;
rxb = R.xb.*G.b;
rxc = R.xc.*G.c;
rya = R.ya.*G.a;
ryb = R.yb.*G.b;
ryc = R.yc.*G.c;
rza = R.za.*G.a;
rzb = R.zb.*G.b;
rzc = R.zc.*G.c;

figure(1), clf
subplot(1,2,1)
plot(1:td1,rxa,'ro',1:td1,rxb,'rs',1:td1,rxc,'rx',...
    1:td1,rya,'go',1:td1,ryb,'gs',1:td1,ryc,'gx',...
    1:td1,rza,'bo',1:td1,rzb,'bs',1:td1,rzc,'bx')

ylabel('g ramps')

title(['td1=' num2str(td1)])

Gnorm.x = sv.bT.dir.x;
Gnorm.y = sv.bT.dir.y;
Gnorm.z = sv.bT.dir.z;

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
plot(X,Y,'k.')
hold on
plot(latitude.X,latitude.Y,'b-')
plot(longitude.X,longitude.Y,'b-')
axis tight equal off

%save gradient ramps

sv.bT.ramp = struct('xa',rxa,'ya',rya,'za',rza,'xb',rxb,'yb',ryb,'zb',rzb,'xc',rxc,'yc',ryc,'zc',rzc);

fpath = [DataDir '/' ExpNam '/' num2str(expno)];
eval(['save ' fpath '/SamplingVector sv'])

fname = 'rxa';
fidxa = fopen([fpath '/' fname],'w');
fname = 'rya';
fidya = fopen([fpath '/' fname],'w');
fname = 'rza';
fidza = fopen([fpath '/' fname],'w');
fname = 'rxb';
fidxb = fopen([fpath '/' fname],'w');
fname = 'ryb';
fidyb = fopen([fpath '/' fname],'w');
fname = 'rzb';
fidzb = fopen([fpath '/' fname],'w');
fname = 'rxc';
fidxc = fopen([fpath '/' fname],'w');
fname = 'ryc';
fidyc = fopen([fpath '/' fname],'w');
fname = 'rzc';
fidzc = fopen([fpath '/' fname],'w');

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
    fprintf(fidxa,'%s\n',text.header{nline}{1});
    fprintf(fidya,'%s\n',text.header{nline}{1});
    fprintf(fidza,'%s\n',text.header{nline}{1});
    fprintf(fidxb,'%s\n',text.header{nline}{1});
    fprintf(fidyb,'%s\n',text.header{nline}{1});
    fprintf(fidzb,'%s\n',text.header{nline}{1});
    fprintf(fidxc,'%s\n',text.header{nline}{1});
    fprintf(fidyc,'%s\n',text.header{nline}{1});
    fprintf(fidzc,'%s\n',text.header{nline}{1});
end


fprintf(fidxa,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidxa,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidya,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidya,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidza,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidza,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidxb,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidxb,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidyb,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidyb,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidzb,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidzb,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidxc,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidxc,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidyc,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidyc,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidzc,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidzc,'%s\n',['##XYDATA= (X++(Y..Y))']);

for nG = 1:length(rxa)
    fprintf(fidxa,'%f\n',rxa(nG));
    fprintf(fidya,'%f\n',rya(nG));
    fprintf(fidza,'%f\n',rza(nG));
    fprintf(fidxb,'%f\n',rxb(nG));
    fprintf(fidyb,'%f\n',ryb(nG));
    fprintf(fidzb,'%f\n',rzb(nG));
    fprintf(fidxc,'%f\n',rxc(nG));
    fprintf(fidyc,'%f\n',ryc(nG));
    fprintf(fidzc,'%f\n',rzc(nG));
end


fprintf(fidxa,'%s\n',['##END']);
fprintf(fidya,'%s\n',['##END']);
fprintf(fidza,'%s\n',['##END']);
fprintf(fidxb,'%s\n',['##END']);
fprintf(fidyb,'%s\n',['##END']);
fprintf(fidzb,'%s\n',['##END']);
fprintf(fidxc,'%s\n',['##END']);
fprintf(fidyc,'%s\n',['##END']);
fprintf(fidzc,'%s\n',['##END']);

fclose(fidxa);
fclose(fidya);
fclose(fidza);
fclose(fidxb);
fclose(fidyb);
fclose(fidzb);
fclose(fidxc);
fclose(fidyc);
fclose(fidzc);

ListDir = '/opt/topspin2/exp/stan/nmr/lists/f1';

fpath = [ListDir];

fname = 'DT_fq';
fid = fopen([fpath '/' fname],'w');

fprintf(fid,'%s\n','sfo hz');
for n = 1:td1
    fprintf(fid,'%1.2f\n',sv.fq(n));
end

fclose(fid);


ListDir = '/opt/topspin2/exp/stan/nmr/lists/vd';

fpath = [ListDir];

fname = 'DT_vdT1';
fid1 = fopen([fpath '/' fname],'w');
fname = 'DT_vdT2';
fid2 = fopen([fpath '/' fname],'w');

for n = 1:td1
    fprintf(fid1,'%1.8f\n',sv.vdT1(n));
    fprintf(fid2,'%1.8f\n',sv.vdT2(n));
end

fclose(fid1);
fclose(fid2);

cd(wd)