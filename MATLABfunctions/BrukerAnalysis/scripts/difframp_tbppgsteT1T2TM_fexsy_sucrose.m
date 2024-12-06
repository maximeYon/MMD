clear all

ListDir = '/opt/topspin4.0.7/exp/stan/nmr/lists/vd';
gpDir = '/opt/topspin4.0.7/exp/stan/nmr/lists/gp/user';
% DataDir = '/opt/nmrdata/daniel/';
% ExpNam = '20210301_sucrose'; expno = 4;

G_filt.min = .05;
G_filt.max = .5;

G_det.min = .05;
G_det.max = 1;
G_det.N = 16;

TM.min = 20e-3;
TM.max = 2000e-3;
TM.N = 7;

%----------------------------------------------
G_filt.N = 2; 
G_filt.theta = pi/2; G_filt.phi = 1.1*pi/4;
G_mix.min = .05;
G_mix.max = .05;
G_mix.N = 1;
G_mix.theta = .1*pi/2; G_mix.phi = .1*pi/4;
G_det.theta = pi/2; G_det.phi = .9*pi/4;
TE.min = 2e-3; TE.max = 2e-3; TE.N = 1;
TR.min = 5000e-3; TR.max = 5000e-3; TR.N = 1;

ga = linspace(G_filt.min,G_filt.max,G_filt.N)';
gb = linspace(G_mix.min,G_mix.max,G_mix.N)';
gc = linspace(G_det.min,G_det.max,G_det.N)';
tm = linspace(TM.min,TM.max,TM.N)'; tm = [TM.min; tm];
te = linspace(TE.min,TE.max,TE.N)';
tr = flipud(linspace(TR.min,TR.max,TR.N)');

[gc,tm] = ndgrid(gc,tm);
ga = G_filt.max*ones(G_det.N,TM.N+1);
ga(:,1) = G_filt.min;
gb = G_mix.max*ones(G_det.N,TM.N+1);
te = TE.max*ones(G_det.N,TM.N+1);
tr = TR.max*ones(G_det.N,TM.N+1);

td1 = numel(ga);
ga = reshape(ga,[td1 1]);
gb = reshape(gb,[td1 1]);
gc = reshape(gc,[td1 1]);
tm = reshape(tm,[td1 1]);
te = reshape(te,[td1 1]);
tr = reshape(tr,[td1 1]);

gax = ga*sin(G_filt.theta)*cos(G_filt.phi);
gay = ga*sin(G_filt.theta)*sin(G_filt.phi);
gaz = ga*cos(G_filt.theta);
gbx = gb*sin(G_mix.theta)*cos(G_mix.phi);
gby = gb*sin(G_mix.theta)*sin(G_mix.phi);
gbz = gb*cos(G_mix.theta);
gcx = gc*sin(G_det.theta)*cos(G_det.phi);
gcy = gc*sin(G_det.theta)*sin(G_det.phi);
gcz = gc*cos(G_det.theta);

vdT1 = tr;
vdT2 = te;
vdTM = tm;

figure(2), clf
subplot(2,2,1)
plot(1:td1,vdT1,'o')
subplot(2,2,2)
plot(1:td1,vdT2,'o')
subplot(2,2,3)
plot(1:td1,vdTM,'o')
subplot(2,2,4)
plot(1:td1,gax,'rx',1:td1,gay,'gx',1:td1,gaz,'bx',...
    1:td1,gbx,'rs',1:td1,gby,'gs',1:td1,gbz,'bs',...
    1:td1,gcx,'ro',1:td1,gcy,'go',1:td1,gcz,'bo')
title(num2str(td1))

wd = cd;


fpath = [ListDir];

fname = 'DT_vdT1';
fid1 = fopen([fpath '/' fname],'w');
fname = 'DT_vdT2';
fid2 = fopen([fpath '/' fname],'w');
fname = 'DT_vdTM';
fid3 = fopen([fpath '/' fname],'w');

for n = 1:td1
    fprintf(fid1,'%1.8f\n',vdT1(n));
    fprintf(fid2,'%1.8f\n',vdT2(n));
    fprintf(fid3,'%1.8f\n',vdTM(n));
end

fclose(fid1);
fclose(fid2);
fclose(fid3);

%fpath = [DataDir '/' ExpNam '/' num2str(expno)];
fpath = gpDir;

fname = 'gax';
fidax = fopen([fpath '/' fname],'w');
fname = 'gay';
fiday = fopen([fpath '/' fname],'w');
fname = 'gaz';
fidaz = fopen([fpath '/' fname],'w');
fname = 'gbx';
fidbx = fopen([fpath '/' fname],'w');
fname = 'gby';
fidby = fopen([fpath '/' fname],'w');
fname = 'gbz';
fidbz = fopen([fpath '/' fname],'w');
fname = 'gcx';
fidcx = fopen([fpath '/' fname],'w');
fname = 'gcy';
fidcy = fopen([fpath '/' fname],'w');
fname = 'gcz';
fidcz = fopen([fpath '/' fname],'w');

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

for nG = 1:length(gax)
    fprintf(fidax,'%f\n',gax(nG));
    fprintf(fiday,'%f\n',gay(nG));
    fprintf(fidaz,'%f\n',gaz(nG));
    fprintf(fidbx,'%f\n',gbx(nG));
    fprintf(fidby,'%f\n',gby(nG));
    fprintf(fidbz,'%f\n',gbz(nG));
    fprintf(fidcx,'%f\n',gcx(nG));
    fprintf(fidcy,'%f\n',gcy(nG));
    fprintf(fidcz,'%f\n',gcz(nG));
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
fclose(fidbx);
fclose(fidby);
fclose(fidbz);
fclose(fidcx);
fclose(fidcy);
fclose(fidcz);

cd(wd)