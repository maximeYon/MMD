clear all

wd = cd;
DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
ExpNam = 'pgse4_test'; expno = 16;

NG = 8;
%Grel = linspace(0,1,NG+1); Grel = Grel(2:NG+1);
Grel = linspace(0.01,1,NG);
%Grel = logspace(-2,0,NG);

theta1 = acos(1/sqrt(3)); phi1 = pi/4;
Gnorm1.x = [sin(theta1).*cos(phi1)];
Gnorm1.y = [sin(theta1).*sin(phi1)];
Gnorm1.z = [cos(theta1)];
theta2 = acos(1/sqrt(3)); phi2 = pi/4+pi;
Gnorm2.x = [sin(theta2).*cos(phi2)];
Gnorm2.y = [sin(theta2).*sin(phi2)];
Gnorm2.z = [cos(theta2)];
theta3 = pi-acos(1/sqrt(3)); phi3 = -pi/4;
Gnorm3.x = [sin(theta3).*cos(phi3)];
Gnorm3.y = [sin(theta3).*sin(phi3)];
Gnorm3.z = [cos(theta3)];
theta4 = pi-acos(1/sqrt(3)); phi4 = pi-pi/4;
Gnorm4.x = [sin(theta4).*cos(phi4)];
Gnorm4.y = [sin(theta4).*sin(phi4)];
Gnorm4.z = [cos(theta4)];

figure(1), clf
plot3([0 Gnorm1.x],[0 Gnorm1.y],[0 Gnorm1.z],'-o',...)
[0 Gnorm2.x],[0 Gnorm2.y],[0 Gnorm2.z],'-o',...
[0 Gnorm3.x],[0 Gnorm3.y],[0 Gnorm3.z],'-o',...
[0 Gnorm4.x],[0 Gnorm4.y],[0 Gnorm4.z],'-o')
axis([-1 1 -1 1 -1 1])
axis square

G1.x = [];
G1.y = [];
G1.z = [];
G2.x = [];
G2.y = [];
G2.z = [];
G3.x = [];
G3.y = [];
G3.z = [];
G4.x = [];
G4.y = [];
G4.z = [];

for nG = 1:NG
    G1.x = [G1.x; Grel(nG)*Gnorm1.x];
    G1.y = [G1.y; Grel(nG)*Gnorm1.y];
    G1.z = [G1.z; Grel(nG)*Gnorm1.z];
    G2.x = [G2.x; Grel(nG)*Gnorm2.x];
    G2.y = [G2.y; Grel(nG)*Gnorm2.y];
    G2.z = [G2.z; Grel(nG)*Gnorm2.z];
    G3.x = [G3.x; Grel(nG)*Gnorm3.x];
    G3.y = [G3.y; Grel(nG)*Gnorm3.y];
    G3.z = [G3.z; Grel(nG)*Gnorm3.z];
    G4.x = [G4.x; Grel(nG)*Gnorm4.x];
    G4.y = [G4.y; Grel(nG)*Gnorm4.y];
    G4.z = [G4.z; Grel(nG)*Gnorm4.z];
end

td1 = NG^4;

[G1.x,G2.x,G3.x,G4.x] = ndgrid(G1.x,G2.x,G3.x,G4.x);
[G1.y,G2.y,G3.y,G4.y] = ndgrid(G1.y,G2.y,G3.y,G4.y);
[G1.z,G2.z,G3.z,G4.z] = ndgrid(G1.z,G2.z,G3.z,G4.z);
G1.x = reshape(G1.x,td1,1);
G2.x = reshape(G2.x,td1,1);
G3.x = reshape(G3.x,td1,1);
G4.x = reshape(G4.x,td1,1);
G1.y = reshape(G1.y,td1,1);
G2.y = reshape(G2.y,td1,1);
G3.y = reshape(G3.y,td1,1);
G4.y = reshape(G4.y,td1,1);
G1.z = reshape(G1.z,td1,1);
G2.z = reshape(G2.z,td1,1);
G3.z = reshape(G3.z,td1,1);
G4.z = reshape(G4.z,td1,1);

figure(1), clf
plot(1:td1,G1.x,'o',1:td1,G1.y,'o',1:td1,G1.z,'o')
hold on
plot(1:td1,G2.x,'s',1:td1,G2.y,'s',1:td1,G2.z,'s')
plot(1:td1,G3.x,'v',1:td1,G3.y,'v',1:td1,G3.z,'v')
plot(1:td1,G4.x,'x',1:td1,G4.y,'x',1:td1,G4.z,'x')

fpath = [DataDir '/' ExpNam '/' num2str(expno)];

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
fname = 'rx2';
fidx2 = fopen([fpath '/' fname],'w');
fname = 'rx2.txt';
fidx2txt = fopen([fpath '/' fname],'w');
fname = 'ry2';
fidy2 = fopen([fpath '/' fname],'w');
fname = 'ry2.txt';
fidy2txt = fopen([fpath '/' fname],'w');
fname = 'rz2';
fidz2 = fopen([fpath '/' fname],'w');
fname = 'rz2.txt';
fidz2txt = fopen([fpath '/' fname],'w');
fname = 'rx3';
fidx3 = fopen([fpath '/' fname],'w');
fname = 'rx3.txt';
fidx3txt = fopen([fpath '/' fname],'w');
fname = 'ry3';
fidy3 = fopen([fpath '/' fname],'w');
fname = 'ry3.txt';
fidy3txt = fopen([fpath '/' fname],'w');
fname = 'rz3';
fidz3 = fopen([fpath '/' fname],'w');
fname = 'rz3.txt';
fidz3txt = fopen([fpath '/' fname],'w');
fname = 'rx4';
fidx4 = fopen([fpath '/' fname],'w');
fname = 'rx4.txt';
fidx4txt = fopen([fpath '/' fname],'w');
fname = 'ry4';
fidy4 = fopen([fpath '/' fname],'w');
fname = 'ry4.txt';
fidy4txt = fopen([fpath '/' fname],'w');
fname = 'rz4';
fidz4 = fopen([fpath '/' fname],'w');
fname = 'rz4.txt';
fidz4txt = fopen([fpath '/' fname],'w');

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
    fprintf(fidx1,'%s\n',text.header{nline}{1});
    fprintf(fidy1,'%s\n',text.header{nline}{1});
    fprintf(fidz1,'%s\n',text.header{nline}{1});
    fprintf(fidx2,'%s\n',text.header{nline}{1});
    fprintf(fidy2,'%s\n',text.header{nline}{1});
    fprintf(fidz2,'%s\n',text.header{nline}{1});
    fprintf(fidx3,'%s\n',text.header{nline}{1});
    fprintf(fidy3,'%s\n',text.header{nline}{1});
    fprintf(fidz3,'%s\n',text.header{nline}{1});
    fprintf(fidx4,'%s\n',text.header{nline}{1});
    fprintf(fidy4,'%s\n',text.header{nline}{1});
    fprintf(fidz4,'%s\n',text.header{nline}{1});
end


fprintf(fidx1,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidx1,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy1,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidy1,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz1,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidz1,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidx2,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidx2,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy2,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidy2,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz2,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidz2,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidx3,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidx3,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy3,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidy3,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz3,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidz3,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidx4,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidx4,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy4,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidy4,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz4,'%s\n',['##NPOINTS=' num2str(td1)]);
fprintf(fidz4,'%s\n',['##XYDATA= (X++(Y..Y))']);

for nG = 1:td1
    fprintf(fidx1,'%f\n',G1.x(nG));
    fprintf(fidx1txt,'%f\n',G1.x(nG));
    fprintf(fidy1,'%f\n',G1.y(nG));
    fprintf(fidy1txt,'%f\n',G1.y(nG));
    fprintf(fidz1,'%f\n',G1.z(nG));
    fprintf(fidz1txt,'%f\n',G1.z(nG));
    fprintf(fidx2,'%f\n',G2.x(nG));
    fprintf(fidx2txt,'%f\n',G2.x(nG));
    fprintf(fidy2,'%f\n',G2.y(nG));
    fprintf(fidy2txt,'%f\n',G2.y(nG));
    fprintf(fidz2,'%f\n',G2.z(nG));
    fprintf(fidz2txt,'%f\n',G2.z(nG));
    fprintf(fidx3,'%f\n',G3.x(nG));
    fprintf(fidx3txt,'%f\n',G3.x(nG));
    fprintf(fidy3,'%f\n',G3.y(nG));
    fprintf(fidy3txt,'%f\n',G3.y(nG));
    fprintf(fidz3,'%f\n',G3.z(nG));
    fprintf(fidz3txt,'%f\n',G3.z(nG));
    fprintf(fidx4,'%f\n',G4.x(nG));
    fprintf(fidx4txt,'%f\n',G4.x(nG));
    fprintf(fidy4,'%f\n',G4.y(nG));
    fprintf(fidy4txt,'%f\n',G4.y(nG));
    fprintf(fidz4,'%f\n',G4.z(nG));
    fprintf(fidz4txt,'%f\n',G4.z(nG));
end

fprintf(fidx1,'%s\n',['##END']);
fprintf(fidy1,'%s\n',['##END']);
fprintf(fidz1,'%s\n',['##END']);
fprintf(fidx2,'%s\n',['##END']);
fprintf(fidy2,'%s\n',['##END']);
fprintf(fidz2,'%s\n',['##END']);
fprintf(fidx3,'%s\n',['##END']);
fprintf(fidy3,'%s\n',['##END']);
fprintf(fidz3,'%s\n',['##END']);
fprintf(fidx4,'%s\n',['##END']);
fprintf(fidy4,'%s\n',['##END']);
fprintf(fidz4,'%s\n',['##END']);

fclose(fidx1);
fclose(fidy1);
fclose(fidz1);
fclose(fidx1txt);
fclose(fidy1txt);
fclose(fidz1txt);
fclose(fidx2);
fclose(fidy2);
fclose(fidz2);
fclose(fidx2txt);
fclose(fidy2txt);
fclose(fidz2txt);
fclose(fidx3);
fclose(fidy3);
fclose(fidz3);
fclose(fidx3txt);
fclose(fidy3txt);
fclose(fidz3txt);
fclose(fidx4);
fclose(fidy4);
fclose(fidz4);
fclose(fidx4txt);
fclose(fidy4txt);
fclose(fidz4txt);
