clear all

wd = cd;
DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%ExpNam = 'imse3d_test'; expno = 28;
%ExpNam = 'AOToct3'; expno = 7;
%ExpNam = 'AOToct4'; expno = 11;
%ExpNam = 'Tissue_qVAS'; expno = 12;
ExpNam = 'AOToct5'; expno = 11;
ExpNam = 'MIC5_2H_setup'; expno = 21;

Grel = [.5 1]; NG = 2;

Gnorm.x = [1 0 0 1/sqrt(2) 1/sqrt(2) 0         1/sqrt(3)]';
Gnorm.y = [0 1 0 1/sqrt(2) 0         1/sqrt(2) 1/sqrt(3)]';
Gnorm.z = [0 0 1 0         1/sqrt(2) 1/sqrt(2) 1/sqrt(3)]';

G.x = [];
G.y = [];
G.z = [];

for nG = 1:length(Grel)
    G.x = [G.x; Grel(nG)*Gnorm.x];
    G.y = [G.y; Grel(nG)*Gnorm.y];
    G.z = [G.z; Grel(nG)*Gnorm.z];
end

G.x = [0.25; 0; G.x];
G.y = [0; 0.25; G.y];
G.z = [0; 0; G.z];

td1 = length(G.x);
figure(1), clf
plot(1:td1,G.x,1:td1,G.y,1:td1,G.z)

fpath = [DataDir '/' ExpNam '/' num2str(expno)];

fname = 'difframp_x';
fidx = fopen([fpath '/' fname],'w');
fname = 'difframp_x.txt';
fidxtxt = fopen([fpath '/' fname],'w');
fname = 'difframp_y';
fidy = fopen([fpath '/' fname],'w');
fname = 'difframp_y.txt';
fidytxt = fopen([fpath '/' fname],'w');
fname = 'difframp_z';
fidz = fopen([fpath '/' fname],'w');
fname = 'difframp_z.txt';
fidztxt = fopen([fpath '/' fname],'w');

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
    fprintf(fidx,'%s\n',text.header{nline}{1});
    fprintf(fidy,'%s\n',text.header{nline}{1});
    fprintf(fidz,'%s\n',text.header{nline}{1});
end


fprintf(fidx,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidx,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidy,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidz,'%s\n',['##XYDATA= (X++(Y..Y))']);

for nG = 1:length(G.x)
    fprintf(fidx,'%f\n',G.x(nG));
    fprintf(fidxtxt,'%f\n',G.x(nG));
    fprintf(fidy,'%f\n',G.y(nG));
    fprintf(fidytxt,'%f\n',G.y(nG));
    fprintf(fidz,'%f\n',G.z(nG));
    fprintf(fidztxt,'%f\n',G.z(nG));
end


fprintf(fidx,'%s\n',['##END']);
fprintf(fidy,'%s\n',['##END']);
fprintf(fidz,'%s\n',['##END']);

fclose(fidx);
fclose(fidy);
fclose(fidz);
fclose(fidxtxt);
fclose(fidytxt);
fclose(fidztxt);