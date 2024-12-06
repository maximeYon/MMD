clear all

wd = cd;
fpath = wd;

DataDir = '/opt/topspin2/data/DT/nmr';
% %DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%ExpNam = 'qMASopt2'; expno = 4;
%ExpNam = 'uFA_RARE'; expno = 101;
ExpNam = 'DiffVACSY'; expno = 55;
fpath = [DataDir '/' ExpNam '/' num2str(expno)];

dt = 1;

fid = fopen('x.txt');
gx = fscanf(fid,'%f');
fid = fopen('y.txt');
gy = fscanf(fid,'%f');
fid = fopen('z.txt');
gz = fscanf(fid,'%f');
fclose(fid);
                
qx = cumsum(gx*dt);
qy = cumsum(gy*dt);
qz = cumsum(gz*dt);
qr = sqrt(qx.^2 + qy.^2 + qz.^2);
gr = gradient(qr)/dt;
maxg = max(abs([gx; gy; gz; gr]));

G.x = gx/maxg;
G.y = gy/maxg;
G.z = gz/maxg;
G.r = gr/maxg;

Nt = length(gx);
t = 1:Nt;

lw = 5;
figure(1), clf
axh1 = axes('position',[0 0 1 .5])
ph1 = plot(t,G.x,t,G.y,t,G.z,t,G.r);
axis tight
axis off
set(ph1,'LineStyle','-',{'Color'},{[1 0 0];[0 .6 0];[0 0 1];[0 0 0]})
axh2 = axes('position',[0 .5 1 .5])
ph2 = plot(t,qx,t,qy,t,qz,t,qr);
axis tight
axis off
set(ph2,'LineStyle','-',{'Color'},{[1 0 0];[0 .6 0];[0 0 1];[0 0 0]})


fname = ['qMASx'];
fidx = fopen([fpath '/' fname],'w');
fname = 'qMASx.txt';
fidxtxt = fopen([fpath '/' fname],'w');
fname = ['qMASy'];
fidy = fopen([fpath '/' fname],'w');
fname = 'qMASy.txt';
fidytxt = fopen([fpath '/' fname],'w');
fname = ['qMASz'];
fidz = fopen([fpath '/' fname],'w');
fname = 'qMASz.txt';
fidztxt = fopen([fpath '/' fname],'w');
fname = ['qMASr'];
fidr = fopen([fpath '/' fname],'w');
fname = 'qMASr.txt';
fidrtxt = fopen([fpath '/' fname],'w');

text.header = {
{'##TITLE= Optimized q-MAS'};
{'##JCAMP-DX= 5.00 Bruker JCAMP library'}
{'##DATA TYPE= Shape Data'}
{'##ORIGIN= Bruker Analytik GmbH'}
{'##OWNER= <nmrsu>'}
{'##DATE= xx'}
{'##TIME= xx'}
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
    fprintf(fidr,'%s\n',text.header{nline}{1});
end


fprintf(fidx,'%s\n',['##NPOINTS=' num2str(Nt)]);
fprintf(fidx,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy,'%s\n',['##NPOINTS=' num2str(Nt)]);
fprintf(fidy,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz,'%s\n',['##NPOINTS=' num2str(Nt)]);
fprintf(fidz,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidr,'%s\n',['##NPOINTS=' num2str(Nt)]);
fprintf(fidr,'%s\n',['##XYDATA= (X++(Y..Y))']);

for nG = 1:length(G.x)
    fprintf(fidx,'%f\n',G.x(nG));
    fprintf(fidxtxt,'%f\n',G.x(nG));
    fprintf(fidy,'%f\n',G.y(nG));
    fprintf(fidytxt,'%f\n',G.y(nG));
    fprintf(fidz,'%f\n',G.z(nG));
    fprintf(fidztxt,'%f\n',G.z(nG));
    fprintf(fidr,'%f\n',G.r(nG));
    fprintf(fidrtxt,'%f\n',G.r(nG));
end


fprintf(fidx,'%s\n',['##END']);
fprintf(fidy,'%s\n',['##END']);
fprintf(fidz,'%s\n',['##END']);
fprintf(fidr,'%s\n',['##END']);

fclose(fidx);
fclose(fidy);
fclose(fidz);
fclose(fidr);
fclose(fidxtxt);
fclose(fidytxt);
fclose(fidztxt);
fclose(fidrtxt);
