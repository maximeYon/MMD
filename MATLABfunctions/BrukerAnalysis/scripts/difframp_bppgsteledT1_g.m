clear all

ListDir = '/opt/topspin4.0.7/exp/stan/nmr/lists/vd';
gpDir = '/opt/topspin4.0.7/exp/stan/nmr/lists/gp/user';
list_prefix = 'DT_';

G_det.min = 0.01; % pgse gradient 
G_det.max = 1;
G_det.N = 64;
G_det.theta = pi/2; G_det.phi = pi/4;

TM = 10e-3; % consstant TM
TR = 1000e-3; % constant TR

%gc = linspace(G_det.min,G_det.max,G_det.N)';
gc = logspace(log10(G_det.min),log10(G_det.max),G_det.N)';

td1 = numel(gc);
tm = TM*ones([td1 1]);
tr = TR*ones([td1 1]);

gcx = gc*sin(G_det.theta)*cos(G_det.phi);
gcy = gc*sin(G_det.theta)*sin(G_det.phi);
gcz = gc*cos(G_det.theta);

vdT1 = tr;
vdTM = tm;

%% write the lists
td1 = numel(gcx);

figure(2), clf
subplot(3,1,1)
plot(1:td1,vdT1,'o')
title(num2str(td1))
ylabel('vdT1')
subplot(3,1,2)
semilogy(1:td1,vdTM,'o')
ylabel('vdTM')
subplot(3,1,3)
plot(1:td1,gcx,'ro',1:td1,gcy,'go',1:td1,gcz,'bo')
ylabel('g')
xlabel('td1')

wd = cd;

fpath = [ListDir];

fname = [list_prefix 'vdT1'];
fid1 = fopen([fpath '/' fname],'w');
fname = [list_prefix 'vdTM'];
fid2 = fopen([fpath '/' fname],'w');

for n = 1:td1
    fprintf(fid1,'%1.8f\n',vdT1(n));
    fprintf(fid2,'%1.8f\n',vdTM(n));
end

fclose(fid1);
fclose(fid2);

fpath = gpDir;

fname = [list_prefix 'gcx'];
fidcx = fopen([fpath '/' fname],'w');
fname = [list_prefix 'gcy'];
fidcy = fopen([fpath '/' fname],'w');
fname = [list_prefix 'gcz'];
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
    fprintf(fidcx,'%s\n',text.header{nline}{1});
    fprintf(fidcy,'%s\n',text.header{nline}{1});
    fprintf(fidcz,'%s\n',text.header{nline}{1});
end


fprintf(fidcx,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidcx,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidcy,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidcy,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidcz,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidcz,'%s\n',['##XYDATA= (X++(Y..Y))']);

for nG = 1:length(gcx)
    fprintf(fidcx,'%f\n',gcx(nG));
    fprintf(fidcy,'%f\n',gcy(nG));
    fprintf(fidcz,'%f\n',gcz(nG));
end


fprintf(fidcx,'%s\n',['##END']);
fprintf(fidcy,'%s\n',['##END']);
fprintf(fidcz,'%s\n',['##END']);

fclose(fidcx);
fclose(fidcy);
fclose(fidcz);

cd(wd)