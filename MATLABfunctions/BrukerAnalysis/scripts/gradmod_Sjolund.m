clear all

wd = cd;
fpath = wd;

DataDir = '/opt/topspin2/data/DT/nmr';
% %DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%ExpNam = 'qMASopt2'; expno = 4;
%ExpNam = 'C14E5_3'; expno = 14;
ExpNam = 'SjolundOpt'; expno = 8;
fpath = [DataDir '/' ExpNam '/' num2str(expno)];

dt = 1;

fid = fopen('SjolundOpt1x.txt');
gx = fscanf(fid,'%f');
fid = fopen('SjolundOpt1y.txt');
gy = fscanf(fid,'%f');
fid = fopen('SjolundOpt1z.txt');
gz = fscanf(fid,'%f');
fclose(fid);

t = linspace(0,1,1000)';
t1 = linspace(0,1,length(gx))';

gx = interp1(t1,gx,t);
gy = interp1(t1,gy,t);
gz = interp1(t1,gz,t);

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

bmat.xx = sum(qx.*qx);
bmat.yy = sum(qy.*qy);
bmat.zz = sum(qz.*qz);
bmat.xy = sum(qx.*qy);
bmat.xz = sum(qx.*qz);
bmat.yz = sum(qy.*qz);

Nt = length(gx);
t = 1:Nt;

lw = 5;
figure(1), clf
axh1 = axes('position',[0.1 0.1 .85 .4])
ph1 = plot(t/Nt,G.x,t/Nt,G.y,t/Nt,G.z,t/Nt,G.r);
axis tight
%axis off
set(ph1,'LineStyle','-',{'Color'},{[1 0 0];[0 .6 0];[0 0 1];[0 0 0]})
xlabel('t/\tau')
ylabel('G / Gmax')

axh2 = axes('position',[0.1 .55 .85 .4])
ph2 = plot(t/Nt,qx/max(qr),t/Nt,qy/max(qr),t/Nt,qz/max(qr),t/Nt,qr/max(qr));
axis tight
%axis off
set(ph2,'LineStyle','-',{'Color'},{[1 0 0];[0 .6 0];[0 0 1];[0 0 0]})
ylabel('q / qmax')


d.t = t/Nt;
lw = 30;
figure(2), clf
axh1 = axes('position',[0 0 1 1]);
ph1 = plot(d.t,G.x,d.t,G.y,d.t,G.z);
axis tight
axis off
set(ph1,'LineStyle','-',{'Color'},{[1 0 0];[0 .6 0];[.1 .1 1]})
set(ph1(1),'LineWidth',lw)
set(ph1(2),'LineWidth',.75*lw)
set(ph1(3),'LineWidth',.5*lw)
set(axh1,'YLim',1.2*[-1 1],'XLim',[-.1 1.1])

set(gcf, 'PaperPosition', [0 0 11 4],'PaperSize', [11 4]); 
eval(['print Sjolund_GxGyGz.pdf -dpdf'])

figure(4), clf
axh1 = axes('position',[0 0 1 1]);
ph1 = plot(d.t,G.r);
axis tight
axis off
set(ph1,'LineStyle','-',{'Color'},{[0 0 0]})
set(ph1(1),'LineWidth',lw)
set(axh1,'YLim',1.2*[-1 1],'XLim',[-.1 1.1])

set(gcf, 'PaperPosition', [0 0 11 4],'PaperSize', [11 4]); 
eval(['print Sjolund_Gdir -dpdf'])
return
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
