clear all

wd = cd;
DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
ExpNam = 'MIC5_2H_setup'; expno = 18;

NG = 4;
%Grel = linspace(0,1,NG+1); Grel = Grel(2:NG+1);
Grel = linspace(0.01,1,NG);
%Grel = logspace(-2,0,NG);

load UniformDistSphereRepulsionN15

Ndir = 2*length(theta);

Gnorm.x = zeros(Ndir,1);
Gnorm.y = zeros(Ndir,1);
Gnorm.z = zeros(Ndir,1);
Gnorm.iso = zeros(Ndir,1);

Gnorm.x(2:2:Ndir) = sin(theta).*cos(phi);
Gnorm.y(2:2:Ndir) = sin(theta).*sin(phi);
Gnorm.z(2:2:Ndir) = cos(theta);
Gnorm.iso(1:2:(Ndir-1)) = 1;

% Gnorm.x = [0; sin(theta).*cos(phi)];
% Gnorm.y = [0; sin(theta).*sin(phi)];
% Gnorm.z = [0; cos(theta)];
% Gnorm.iso = [1; zeros(size(theta))];

% Gnorm.x = [0 1 0 0 1/sqrt(2) 1/sqrt(2) 0         1/sqrt(3)]';
% Gnorm.y = [0 0 1 0 1/sqrt(2) 0         1/sqrt(2) 1/sqrt(3)]';
% Gnorm.z = [0 0 0 1 0         1/sqrt(2) 1/sqrt(2) 1/sqrt(3)]';
% Gnorm.iso = [1 0 0 0 0 0 0 0]';


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

td1 = length(G.x);

thresh = 0.01;
%G.x(find(G.x<thresh)) = thresh;
%G.y(find(G.y<thresh)) = thresh;
%G.z(find(G.z<thresh)) = thresh;

figure(1), clf
subplot(1,2,1)
plot(1:td1,G.x,'ro',1:td1,G.y,'go',1:td1,G.z,'bo',1:td1,G.iso,'ko')

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
end


fprintf(fidx,'%s\n',['##NPOINTS=' num2str(NG*length(Gnorm.x))]);
fprintf(fidx,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy,'%s\n',['##NPOINTS=' num2str(NG*length(Gnorm.x))]);
fprintf(fidy,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz,'%s\n',['##NPOINTS=' num2str(NG*length(Gnorm.x))]);
fprintf(fidz,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidiso,'%s\n',['##NPOINTS=' num2str(NG*length(Gnorm.x))]);
fprintf(fidiso,'%s\n',['##XYDATA= (X++(Y..Y))']);

for nG = 1:length(G.x)
    fprintf(fidx,'%f\n',G.x(nG));
    fprintf(fidxtxt,'%f\n',G.x(nG));
    fprintf(fidy,'%f\n',G.y(nG));
    fprintf(fidytxt,'%f\n',G.y(nG));
    fprintf(fidz,'%f\n',G.z(nG));
    fprintf(fidztxt,'%f\n',G.z(nG));
    fprintf(fidiso,'%f\n',G.iso(nG));
    fprintf(fidisotxt,'%f\n',G.iso(nG));
end


fprintf(fidx,'%s\n',['##END']);
fprintf(fidy,'%s\n',['##END']);
fprintf(fidz,'%s\n',['##END']);
fprintf(fidiso,'%s\n',['##END']);

fclose(fidx);
fclose(fidy);
fclose(fidz);
fclose(fidxtxt);
fclose(fidytxt);
fclose(fidztxt);
fclose(fidiso);
fclose(fidisotxt);

td1