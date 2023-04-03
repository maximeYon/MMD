clear all

wd = cd;
DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
ExpNam = 'qMASopt'; expno = 90;
ExpNam = 'qMASdtirare2d_test'; expno = 13;
ExpNam = 'isoanisocorr_qMASopt'; expno = 39;

NG = 16;
%Grel = linspace(0,1,NG+1); Grel = Grel(2:NG+1);
Grel = linspace(0.01,1,NG);
%Grel = logspace(-2,0,NG);

load UniformDistSphereRepulsionN15

Gnorm.x = [0; sin(theta).*cos(phi)];
Gnorm.y = [0; sin(theta).*sin(phi)];
Gnorm.z = [0; cos(theta)];
Gnorm.iso = [1; zeros(size(theta))];

%Gnorm.x = [0 1 0 0 1/sqrt(2) 1/sqrt(2) 0         1/sqrt(3)]';
%Gnorm.y = [0 0 1 0 1/sqrt(2) 0         1/sqrt(2) 1/sqrt(3)]';
%Gnorm.z = [0 0 0 1 0         1/sqrt(2) 1/sqrt(2) 1/sqrt(3)]';
%Gnorm.iso = [1 0 0 0 0 0 0 0]';




[Gnorm.x, Gnorm1.x] = ndgrid(Gnorm.x,Gnorm.x);
Gnorm.x = reshape(Gnorm.x,numel(Gnorm.x),1)
Gnorm1.x = reshape(Gnorm1.x,numel(Gnorm1.x),1)
[Gnorm.y, Gnorm1.y] = ndgrid(Gnorm.y,Gnorm.y);
Gnorm.y = reshape(Gnorm.y,numel(Gnorm.y),1)
Gnorm1.y = reshape(Gnorm1.y,numel(Gnorm1.y),1)
[Gnorm.z, Gnorm1.z] = ndgrid(Gnorm.z,Gnorm.z);
Gnorm.z = reshape(Gnorm.z,numel(Gnorm.z),1)
Gnorm1.z = reshape(Gnorm1.z,numel(Gnorm1.z),1)
[Gnorm.iso, Gnorm1.iso] = ndgrid(Gnorm.iso,Gnorm.iso);
Gnorm.iso = reshape(Gnorm.iso,numel(Gnorm.iso),1)
Gnorm1.iso = reshape(Gnorm1.iso,numel(Gnorm1.iso),1)

NDir2 = numel(Gnorm.x);
figure(1), clf
plot(1:NDir2,Gnorm.x,'o',1:NDir2,Gnorm.y,'o',1:NDir2,Gnorm.z,'o',1:NDir2,Gnorm.iso,'ko')
hold on
plot(1:NDir2,Gnorm1.x,'x',1:NDir2,Gnorm1.y,'x',1:NDir2,Gnorm1.z,'x',1:NDir2,Gnorm1.iso,'kx')

G.x = [];
G.y = [];
G.z = [];
G.iso = [];
G1.x = [];
G1.y = [];
G1.z = [];
G1.iso = [];

td1 = NG*NDir2

for nG = 1:NG
    G.x = [G.x; Grel(nG)*Gnorm.x];
    G.y = [G.y; Grel(nG)*Gnorm.y];
    G.z = [G.z; Grel(nG)*Gnorm.z];
    G.iso = [G.iso; Grel(nG)*Gnorm.iso];
    G1.x = [G1.x; Grel(nG)*Gnorm1.x];
    G1.y = [G1.y; Grel(nG)*Gnorm1.y];
    G1.z = [G1.z; Grel(nG)*Gnorm1.z];
    G1.iso = [G1.iso; Grel(nG)*Gnorm1.iso];
end

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
