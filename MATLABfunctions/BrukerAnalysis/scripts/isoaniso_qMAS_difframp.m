clear all

wd = cd;
DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
ExpNam = 'qMASopt'; expno = 90;
ExpNam = 'qMASdtirare2d_test'; expno = 13;
ExpNam = 'isoanisocorr_qMASopt'; expno = 6;

NG = 16;
Grel = linspace(0,1,NG+1); Grel = Grel(2:NG+1);
%Grel = linspace(0.01,1,NG);
%Grel = logspace(-2,0,NG);

Gnorm.x = [0 1 0 0 1/sqrt(2) 1/sqrt(2) 0         1/sqrt(3)]';
Gnorm.y = [0 0 1 0 1/sqrt(2) 0         1/sqrt(2) 1/sqrt(3)]';
Gnorm.z = [0 0 0 1 0         1/sqrt(2) 1/sqrt(2) 1/sqrt(3)]';
Gnorm.iso = [0 0 0 0 0 0 0 0]';

G.x = [];
G.y = [];
G.z = [];
G.iso = [];

for nG = 1:length(Grel)
    G.x = [G.x; Grel(nG)*Gnorm.x];
    G.y = [G.y; Grel(nG)*Gnorm.y];
    G.z = [G.z; Grel(nG)*Gnorm.z];
    G.iso = [G.iso; Grel(nG)*Gnorm.iso];
end

Gdim1.x = G.x;
Gdim1.y = G.y;
Gdim1.z = G.z;
Gdim1.iso = G.iso;

NG1 = 16;
td1 = NG*8*NG1
G1rel = 1*linspace(0.01,1,NG1);

G1norm.x = zeros(NG*8,1);
G1norm.y = zeros(NG*8,1);
G1norm.z = zeros(NG*8,1);
G1norm.iso = ones(NG*8,1);

G1.x = G1rel(1)*G1norm.x;
G1.y = G1rel(1)*G1norm.y;
G1.z = G1rel(1)*G1norm.z;
G1.iso = G1rel(1)*G1norm.iso;

for nG1 = 2:length(G1rel)
    G1.x = [G1.x; G1rel(nG1)*G1norm.x];
    G1.y = [G1.y; G1rel(nG1)*G1norm.y];
    G1.z = [G1.z; G1rel(nG1)*G1norm.z];
    G1.iso = [G1.iso; G1rel(nG1)*G1norm.iso];
    G.x = [G.x; Gdim1.x];
    G.y = [G.y; Gdim1.y];
    G.z = [G.z; Gdim1.z];
    G.iso = [G.iso; Gdim1.iso];
end

thresh = 0.01;
%G.x(find(G.x<thresh)) = thresh;
%G.y(find(G.y<thresh)) = thresh;
%G.z(find(G.z<thresh)) = thresh;


fpath = [DataDir '/' ExpNam '/' num2str(expno)];

fname = 'rx';
fidx = fopen([fpath '/' fname],'w');
fname = 'rx.txt';
fidxtxt = fopen([fpath '/' fname],'w');
fname = 'ry';
fidy = fopen([fpath '/' fname],'w');
fname = 'ry.txt';
fidytxt = fopen([fpath '/' fname],'w');
fname = 'rz';
fidz = fopen([fpath '/' fname],'w');
fname = 'rz.txt';
fidztxt = fopen([fpath '/' fname],'w');
fname = 'riso';
fidiso = fopen([fpath '/' fname],'w');
fname = 'riso.txt';
fidisotxt = fopen([fpath '/' fname],'w');
fname = 'rx1';
fidx1 = fopen([fpath '/' fname],'w');
fname = 'rx1.txt';
fidx1txt = fopen([fpath '/' fname],'w');
fname = 'ry1';
fidy1 = fopen([fpath '/' fname],'w');
fname = 'ry1.txt';
fidy1txt = fopen([fpath '/' fname],'w');
fname = 'rz1';
fidz1 = fopen([fpath '/' fname],'w');
fname = 'rz1.txt';
fidz1txt = fopen([fpath '/' fname],'w');
fname = 'riso1';
fidiso1 = fopen([fpath '/' fname],'w');
fname = 'riso1.txt';
fidiso1txt = fopen([fpath '/' fname],'w');

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
    fprintf(fidiso,'%s\n',text.header{nline}{1});
    fprintf(fidx1,'%s\n',text.header{nline}{1});
    fprintf(fidy1,'%s\n',text.header{nline}{1});
    fprintf(fidz1,'%s\n',text.header{nline}{1});
    fprintf(fidiso1,'%s\n',text.header{nline}{1});
end


fprintf(fidx,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidx,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidy,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidz,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidiso,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidiso,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidx1,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidx1,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy1,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidy1,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz1,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidz1,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidiso1,'%s\n',['##NPOINTS=' num2str(length(G.x))]);
fprintf(fidiso1,'%s\n',['##XYDATA= (X++(Y..Y))']);

for nG = 1:length(G.x)
    fprintf(fidx,'%f\n',G.x(nG));
    fprintf(fidxtxt,'%f\n',G.x(nG));
    fprintf(fidy,'%f\n',G.y(nG));
    fprintf(fidytxt,'%f\n',G.y(nG));
    fprintf(fidz,'%f\n',G.z(nG));
    fprintf(fidztxt,'%f\n',G.z(nG));
    fprintf(fidiso,'%f\n',G.iso(nG));
    fprintf(fidisotxt,'%f\n',G.iso(nG));
    fprintf(fidx1,'%f\n',G1.x(nG));
    fprintf(fidx1txt,'%f\n',G1.x(nG));
    fprintf(fidy1,'%f\n',G1.y(nG));
    fprintf(fidy1txt,'%f\n',G1.y(nG));
    fprintf(fidz1,'%f\n',G1.z(nG));
    fprintf(fidz1txt,'%f\n',G1.z(nG));
    fprintf(fidiso1,'%f\n',G1.iso(nG));
    fprintf(fidiso1txt,'%f\n',G1.iso(nG));
end


fprintf(fidx,'%s\n',['##END']);
fprintf(fidy,'%s\n',['##END']);
fprintf(fidz,'%s\n',['##END']);
fprintf(fidiso,'%s\n',['##END']);
fprintf(fidx1,'%s\n',['##END']);
fprintf(fidy1,'%s\n',['##END']);
fprintf(fidz1,'%s\n',['##END']);
fprintf(fidiso1,'%s\n',['##END']);

fclose(fidx);
fclose(fidy);
fclose(fidz);
fclose(fidxtxt);
fclose(fidytxt);
fclose(fidztxt);
fclose(fidiso);
fclose(fidisotxt);
fclose(fidx1);
fclose(fidy1);
fclose(fidz1);
fclose(fidx1txt);
fclose(fidy1txt);
fclose(fidz1txt);
fclose(fidiso1);
fclose(fidiso1txt);
