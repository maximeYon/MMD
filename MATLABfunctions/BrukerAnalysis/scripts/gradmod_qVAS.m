clear all

wd = cd;
fpath = wd;

DataDir = '/opt/topspin2/data/DT/nmr';
% %DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%ExpNam = 'DiffVACSY'; expno = 124;
ExpNam = 'C14E5_5'; expno = 11; 

fpath = [DataDir '/' ExpNam '/' num2str(expno)];

noptexport = 21; %number of terms in cos expansion, 21
d.Nt = 1000; %number of time steps

d.tE = 1;

d.q = 1;
d.zetam = pi/2;
d.n = 1;

d.t = d.tE*linspace(0,1,d.Nt)';
d.dt = d.t(2)-d.t(1);

OptDat.Pin{21} = [2.7600e-01
   1.8636e-01
   1.0271e-01
   2.8803e-02
   1.3061e-02
   3.8476e-03
   5.0764e-03
   7.8710e-03
   1.0979e-02
   7.2733e-03
   4.6788e-03
   2.0443e-03
   2.2196e-04
   1.8141e-03
   5.0880e-08
   1.0098e-03
   1.2173e-03
   8.9312e-04
   7.9497e-04
  -2.2231e-04]'; %20, 9.6113

nopt_v = [21];


for nopt = nopt_v
    Pin = OptDat.Pin{nopt};

    psi0 = 0;
    theta = 0;
    phi = 0;
    aF = Pin(1:length(Pin));
    aF = [1-sum(aF) aF];

    F = zeros(size(d.t));

    for na = 1:length(aF)
        F = F + aF(na)*.5*(1-cos(na*2*pi*d.t/d.tE));
    end

    F = F/max(F);

    qt = d.q*F;
    td = sum(F.^2*d.dt);
    psirate = 2*pi*d.n/td*F.^2;
    psi = psi0 + cumsum(psirate*d.dt);

    qx = qt.*sin(d.zetam).*cos(psi);
    qy = qt.*sin(d.zetam).*sin(psi);
    qz = qt.*cos(d.zetam);

    gx = gradient(qx)/d.dt;
    gy = gradient(qy)/d.dt;
    gz = gradient(qz)/d.dt;

    gxold = gx; gzold = gz;
    gx = gxold*cos(theta) + gzold*sin(theta);
    gz = gzold*cos(theta) - gxold*sin(theta);

    gxold = gx; gyold = gy;
    gx = gxold*cos(phi) - gyold*sin(phi);
    gy = gyold*cos(phi) + gxold*sin(phi);
% 
%     gxold = gx; gyold = gy; gzold = gz;
%     gx = gxold; gy = -gzold; gz = gyold;

    gr = sqrt(gx.^2 + gy.^2 + gz.^2);
    gt = gradient(qt)/d.dt;

    qx = cumsum(gx*d.dt);
    qy = cumsum(gy*d.dt);
    qz = cumsum(gz*d.dt);
    qr = sqrt(qx.^2 + qy.^2 + qz.^2);

    b = sum(qr.^2*d.dt);
    maxg = max(abs([gx; gy; gz]));
    C = b/maxg^2;

    OptDat.G.x{nopt} = gx/maxg;
    OptDat.G.y{nopt} = gy/maxg;
    OptDat.G.z{nopt} = gz/maxg;
    OptDat.G.r{nopt} = gt/maxg;
    
    
    OptDat.C{nopt} = C;
    OptDat.energy{nopt} = 1/d.Nt/3*(sum((gx/maxg).^2) + sum((gy/maxg).^2) + sum((gz/maxg).^2));
end


nopt = noptexport;
G.x = OptDat.G.x{nopt};
G.y = OptDat.G.y{nopt};
G.z = OptDat.G.z{nopt};
G.r = OptDat.G.r{nopt};

lw = 5;
figure(1), clf
axh1 = axes('position',[0 0 1 1])
ph1 = plot(d.t,G.x,d.t,G.y,d.t,G.z,d.t,G.r);
axis tight
axis off
set(ph1,'LineStyle','-',{'Color'},{[1 0 0];[0 .6 0];[0 0 1];[0 0 0]})


fname = ['qVASa'];
fidx = fopen([fpath '/' fname],'w');
fname = 'qVASa.txt';
fidxtxt = fopen([fpath '/' fname],'w');
fname = ['qVASb'];
fidy = fopen([fpath '/' fname],'w');
fname = 'qVASb.txt';
fidytxt = fopen([fpath '/' fname],'w');
fname = ['qVASc'];
fidr = fopen([fpath '/' fname],'w');
fname = 'qVASc.txt';
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
    fprintf(fidr,'%s\n',text.header{nline}{1});
end


fprintf(fidx,'%s\n',['##NPOINTS=' num2str(d.Nt)]);
fprintf(fidx,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy,'%s\n',['##NPOINTS=' num2str(d.Nt)]);
fprintf(fidy,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidr,'%s\n',['##NPOINTS=' num2str(d.Nt)]);
fprintf(fidr,'%s\n',['##XYDATA= (X++(Y..Y))']);

for nG = 1:length(G.x)
    fprintf(fidx,'%f\n',G.x(nG));
    fprintf(fidxtxt,'%f\n',G.x(nG));
    fprintf(fidy,'%f\n',G.y(nG));
    fprintf(fidytxt,'%f\n',G.y(nG));
    fprintf(fidr,'%f\n',G.r(nG));
    fprintf(fidrtxt,'%f\n',G.r(nG));
end


fprintf(fidx,'%s\n',['##END']);
fprintf(fidy,'%s\n',['##END']);
fprintf(fidr,'%s\n',['##END']);

fclose(fidx);
fclose(fidy);
fclose(fidr);
fclose(fidxtxt);
fclose(fidytxt);
fclose(fidrtxt);
