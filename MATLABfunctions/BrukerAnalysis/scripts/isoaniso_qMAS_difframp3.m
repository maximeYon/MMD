clear all

wd = cd;
DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
ExpNam = 'qMASopt'; expno = 90;
ExpNam = 'qMASdtirare2d_test'; expno = 13;
ExpNam = 'isoanisocorr_qMASopt'; expno = 41;

NG = 8;
%Grel = linspace(0,1,NG+1); Grel = Grel(2:NG+1);
Grel = linspace(0.01,1,NG);
%Grel = logspace(-2,0,NG);

load UniformDistSphereRepulsionN7

Gnorm.x = [0; sin(theta).*cos(phi)];
Gnorm.y = [0; sin(theta).*sin(phi)];
Gnorm.z = [0; cos(theta)];
Gnorm.iso = [1; zeros(size(theta))];

Gnorm.x = [0 1 0 0 1/sqrt(2) 1/sqrt(2) 0         1/sqrt(3)]';
Gnorm.y = [0 0 1 0 1/sqrt(2) 0         1/sqrt(2) 1/sqrt(3)]';
Gnorm.z = [0 0 0 1 0         1/sqrt(2) 1/sqrt(2) 1/sqrt(3)]';
Gnorm.iso = [1 0 0 0 0 0 0 0]';

% Gnorm.x = [0 1 0 0]';
% Gnorm.y = [0 0 1 0]';
% Gnorm.z = [0 0 0 1]';
% Gnorm.iso = [1 0 0 0]';

G.x = [];
G.y = [];
G.z = [];
G.iso = [];

for nG = 1:NG
    G.x = [G.x; Grel(nG)*Gnorm.x];
    G.y = [G.y; Grel(nG)*Gnorm.y];
    G.z = [G.z; Grel(nG)*Gnorm.z];
    G.iso = [G.iso; Grel(nG)*Gnorm.iso];
end

[G.x, G1.x] = ndgrid(G.x,G.x);
G.x = reshape(G.x,numel(G.x),1)
G1.x = reshape(G1.x,numel(G1.x),1)
[G.y, G1.y] = ndgrid(G.y,G.y);
G.y = reshape(G.y,numel(G.y),1)
G1.y = reshape(G1.y,numel(G1.y),1)
[G.z, G1.z] = ndgrid(G.z,G.z);
G.z = reshape(G.z,numel(G.z),1)
G1.z = reshape(G1.z,numel(G1.z),1)
[G.iso, G1.iso] = ndgrid(G.iso,G.iso);
G.iso = reshape(G.iso,numel(G.iso),1)
G1.iso = reshape(G1.iso,numel(G1.iso),1)

td1 = numel(G.x)


NDir2 = numel(Gnorm.x);
figure(1), clf
plot(1:td1,G.x,'o',1:td1,G.y,'o',1:td1,G.z,'o',1:td1,G.iso,'ko')
hold on
plot(1:td1,G1.x,'x',1:td1,G1.y,'x',1:td1,G1.z,'x',1:td1,G1.iso,'kx')

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
