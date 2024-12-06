clear all

wd = cd;
fpath = wd;

DataDir = '/opt/topspin2/data/DT/nmr';
% %DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%ExpNam = 'qMASopt2'; expno = 4;
%ExpNam = 'uFA_RARE'; expno = 101;
ExpNam = 'DiffVACSY'; expno = 55;
fpath = [DataDir '/' ExpNam '/' num2str(expno)];

noptexport = 21; %number of terms in cos expansion, 1-21
d.Nt = 1000; %number of time steps

d.tE = 1;

d.q = 1;
d.zetam = acos(1/sqrt(3));
d.n = 1;

d.t = d.tE*linspace(0,1,d.Nt)';
d.dt = d.t(2)-d.t(1);

OptDat.Pin{21} = [5.2379e+00
   8.4677e+00
   7.0737e+00
   2.7426e-01
   1.6133e-01
   7.0668e-02
   3.4026e-02
   2.2751e-02
   2.3259e-02
   1.1676e-02
   9.2002e-03
   3.9455e-03
   2.0137e-03
   5.0465e-03
   3.6548e-03
   1.8225e-03
   8.1308e-04
   7.5379e-04
   1.0564e-03
   9.1826e-04
   4.7785e-05
   3.9031e-04
   1.3805e-04]'; %20, 27.8
OptDat.C{21} = 27.77/1000;
OptDat.Pin{20} = [5.2373e+00
   8.4678e+00
   7.0735e+00
   2.7424e-01
   1.6132e-01
   7.0484e-02
   3.3893e-02
   2.2677e-02
   2.3355e-02
   1.1735e-02
   9.1599e-03
   3.9651e-03
   1.9869e-03
   5.1037e-03
   3.7430e-03
   1.9201e-03
   7.6905e-04
   7.4381e-04
   1.0398e-03
   9.7037e-04
   9.5198e-05
   3.9969e-04]'; %19, 27.8
OptDat.C{20} = 27.77/1000;
OptDat.Pin{19} = [5.2353e+00
   8.4694e+00
   7.0703e+00
   2.7469e-01
   1.6193e-01
   7.0767e-02
   3.3266e-02
   2.2606e-02
   2.2951e-02
   1.2125e-02
   9.0037e-03
   3.6811e-03
   2.0978e-03
   4.9030e-03
   3.8877e-03
   1.9409e-03
   6.0821e-04
   6.1521e-04
   1.0807e-03
   9.2408e-04
   1.7800e-04]'; %18, 27.6
OptDat.C{19} = 27.6/1000;
OptDat.Pin{18} = [5.2360e+00
   8.4688e+00
   7.0704e+00
   2.7489e-01
   1.6222e-01
   7.1384e-02
   3.3069e-02
   2.2236e-02
   2.2691e-02
   1.2172e-02
   8.9856e-03
   3.5046e-03
   1.6320e-03
   5.0469e-03
   3.6804e-03
   1.9135e-03
   4.5715e-04
   4.7926e-04
   9.4316e-04
   1.0253e-03]'; %17, 27.6
OptDat.C{18} = 27.56/1000;
OptDat.Pin{17} = [5.2305e+00
   8.4737e+00
   7.0612e+00
   2.7521e-01
   1.6265e-01
   7.1544e-02
   3.2849e-02
   2.2905e-02
   2.1691e-02
   1.1734e-02
   9.2350e-03
   3.5961e-03
   2.1339e-03
   4.8607e-03
   3.5498e-03
   1.6839e-03
   3.9677e-04
   5.1834e-04
   1.2885e-03]'; %16, 27.3
OptDat.C{17} = 27.3/1000;
OptDat.Pin{16} = [5.2324e+00
   8.4712e+00
   7.0660e+00
   2.7587e-01
   1.6173e-01
   7.1735e-02
   3.3287e-02
   2.2757e-02
   2.0473e-02
   1.1643e-02
   8.3681e-03
   3.7837e-03
   1.7832e-03
   4.4313e-03
   2.8218e-03
   1.4962e-03
   4.6735e-04
   8.7021e-04]'; %15, 26.9
OptDat.C{16} = 26.9/1000;
OptDat.Pin{15} = [5.2308e+00
   8.4734e+00
   7.0618e+00
   2.7600e-01
   1.6105e-01
   7.2788e-02
   3.2508e-02
   2.2427e-02
   2.0514e-02
   1.1435e-02
   8.5087e-03
   3.7167e-03
   1.8034e-03
   4.1866e-03
   2.5490e-03
   1.7086e-03
   5.7349e-04]'; %14, 26.8
OptDat.C{15} = 26.8/1000;
OptDat.Pin{14} = [5.2278e+00
   8.4757e+00
   7.0567e+00
   2.7815e-01
   1.6302e-01
   7.2737e-02
   3.1678e-02
   2.2190e-02
   1.9680e-02
   1.0928e-02
   7.7427e-03
   2.7536e-03
   1.1218e-03
   3.7854e-03
   2.1485e-03
   1.1791e-03]'; %13, 26.7
OptDat.C{14} = 26.7/1000;
OptDat.Pin{13} = [5.2406e+00
   8.4649e+00
   7.0796e+00
   2.8175e-01
   1.6126e-01
   7.1031e-02
   3.1534e-02
   2.0129e-02
   2.0357e-02
   1.0878e-02
   7.2446e-03
   1.6984e-03
   6.9696e-04
   4.3223e-03
   2.4292e-03]'; %12, 26.7
OptDat.C{13} = 26.7/1000;
OptDat.Pin{12} = [5.2300e+00
   8.4735e+00
   7.0606e+00
   2.8171e-01
   1.6443e-01
   7.0894e-02
   3.0532e-02
   2.0142e-02
   2.0148e-02
   8.6315e-03
   5.4800e-03
   1.5009e-03
   1.0979e-03
   3.9612e-03]'; %11, 26.4
OptDat.C{12} = 26.4/1000;
OptDat.Pin{11} = [5.2120e+00
   8.4874e+00
   7.0300e+00
   2.8257e-01
   1.6649e-01
   6.5462e-02
   3.1271e-02
   2.4244e-02
   1.7856e-02
   5.1028e-03
   5.3748e-03
   2.5994e-03
   2.3175e-03]'; %10, 26.2
OptDat.C{11} = 26.2/1000;
OptDat.Pin{10} = [5.2063e+00
   8.4913e+00
   7.0211e+00
   2.8402e-01
   1.6689e-01
   6.8836e-02
   2.9798e-02
   2.2664e-02
   1.5068e-02
   5.9755e-03
   4.7838e-03
   1.6639e-03]'; %9, 25.7
OptDat.C{10} = 25.7/1000;
OptDat.Pin{9} = [5.1977e+00
   8.4978e+00
   7.0074e+00
   2.8409e-01
   1.6755e-01
   7.1641e-02
   2.7618e-02
   2.2720e-02
   1.5443e-02
   6.9098e-03
   4.8846e-03]'; %8, 25.4
OptDat.C{9} = 25.42/1000;
OptDat.Pin{8} = [5.1660e+00
   8.5179e+00
   6.9553e+00
   2.8469e-01
   1.7999e-01
   6.6429e-02
   2.5499e-02
   2.1691e-02
   1.4758e-02
   5.9106e-03]'; %7, 25.5
OptDat.C{8} = 25.46/1000;
OptDat.Pin{7} = [5.2035e+00
   8.4924e+00
   7.0185e+00
   2.9540e-01
   1.7274e-01
   5.7517e-02
   2.3512e-02
   1.8395e-02
   1.4822e-02]'; %6, 25.5
OptDat.C{7} = 25.46/1000;
OptDat.Pin{6} = [5.1926e+00
   8.4982e+00
   7.0038e+00
   2.9709e-01
   1.6781e-01
   6.3047e-02
   2.4819e-02
   1.8368e-02]'; %5, 23.2
OptDat.C{6} = 23.2/1000;
OptDat.Pin{5} = [5.2559e+00
   8.4469e+00
   7.1176e+00
   3.0736e-01
   1.4987e-01
   5.2455e-02
   2.5087e-02]'; %4, 22.3
OptDat.C{5} = 22.3/1000;
OptDat.Pin{4} = [5.2109e+00
   8.4752e+00
   7.0555e+00
   3.2086e-01
   1.4212e-01
   3.5750e-02]'; %3, 20.7
OptDat.C{4} = 20.7/1000;
OptDat.Pin{3} = [5.2273e+00
   8.4748e+00
   7.0570e+00
   3.2999e-01
   1.2902e-01]'; %2, 19.4
OptDat.C{3} = 19.4/1000;
OptDat.Pin{2} = [5.2733e+00
   8.4127e+00
   7.1596e+00
   3.3535e-01]'; %1, 12.1
OptDat.C{2} = 12.1/1000;
OptDat.Pin{1} = [5.1419e+00
   8.5325e+00
   6.9152e+00]'; %0, 3.5
OptDat.C{1} = 3.5/1000;


Nopt = 21;


for nopt = 1:Nopt
    Pin = OptDat.Pin{nopt};

    psi0 = Pin(1);
    theta = Pin(2);
    phi = Pin(3);
    aF = Pin(4:length(Pin));
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

    gxold = gx; gyold = gy; gzold = gz;
    gx = gxold; gy = -gzold; gz = gyold;

    gr = sqrt(gx.^2 + gy.^2 + gz.^2);
    gt = gradient(qt)/d.dt;

    qx = cumsum(gx*d.dt);
    qy = cumsum(gy*d.dt);
    qz = cumsum(gz*d.dt);
    qr = sqrt(qx.^2 + qy.^2 + qz.^2);

    b = sum(qr.^2*d.dt);
    maxg = max(abs([gx; gy; gz; gt]));

    OptDat.G.x{nopt} = gx/maxg;
    OptDat.G.y{nopt} = gy/maxg;
    OptDat.G.z{nopt} = gz/maxg;
    OptDat.G.r{nopt} = gt/maxg;
    
    OptDat.energy{nopt} = 1/d.Nt/3*(sum((gx/maxg).^2) + sum((gy/maxg).^2) + sum((gz/maxg).^2));
end

Cvector = zeros(Nopt,1);
for nopt = 1:Nopt
    Cvector(nopt,1) = OptDat.C{nopt};
    Evector(nopt,1) = OptDat.energy{nopt};
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
eval(['print OptiCube_GxGyGz.pdf -dpdf'])

figure(4), clf
axh1 = axes('position',[0 0 1 1]);
ph1 = plot(d.t,G.r);
axis tight
axis off
set(ph1,'LineStyle','-',{'Color'},{[0 0 0]})
set(ph1(1),'LineWidth',lw)
set(axh1,'YLim',1.2*[-1 1],'XLim',[-.1 1.1])

set(gcf, 'PaperPosition', [0 0 11 4],'PaperSize', [11 4]); 
eval(['print OptiCube_Gdir -dpdf'])

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


fprintf(fidx,'%s\n',['##NPOINTS=' num2str(d.Nt)]);
fprintf(fidx,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy,'%s\n',['##NPOINTS=' num2str(d.Nt)]);
fprintf(fidy,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz,'%s\n',['##NPOINTS=' num2str(d.Nt)]);
fprintf(fidz,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidr,'%s\n',['##NPOINTS=' num2str(d.Nt)]);
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
