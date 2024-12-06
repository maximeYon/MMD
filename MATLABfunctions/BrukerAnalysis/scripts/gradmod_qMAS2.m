clear all

wd = cd;
fpath = wd;

DataDir = '/opt/topspin2/data/DT/nmr';
% %DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%ExpNam = 'qMASopt2'; expno = 4;
%ExpNam = 'C14E5_3'; expno = 14;
%ExpNam = 'SjolundOpt'; expno = 9;
ExpNam = 'MIC5_2H_setup'; expno = 18;

fpath = [DataDir '/' ExpNam '/' num2str(expno)];

noptexport = 21; %number of terms in cos expansion, 1-21
d.Nt = 1000; %number of time steps

d.tE = 1;

d.q = 1;
d.zetam = acos(1/sqrt(3));
%d.zetam = pi/2;
d.n = 1;

d.t = d.tE*linspace(0,1,d.Nt)';
d.dt = d.t(2)-d.t(1);

OptDat.Pin{21} = [7.3135e+00
   5.3156e+00
   1.0235e+01
   3.0059e-01
   1.3474e-01
   3.3771e-02
   1.1300e-02
   1.2110e-02
   1.8537e-02
   7.0313e-03
   5.9189e-03
   3.9103e-03
   4.4253e-03
   5.1541e-03
   6.0097e-04
  -2.7835e-05
   1.4314e-03
   1.6354e-03
   8.6584e-04
   4.8698e-04
  -8.6823e-06
   4.7618e-04
   2.2678e-09]'; %20, 7.7366
OptDat.Pin{20}  = [7.3135e+00
   5.3156e+00
   1.0235e+01
   3.0059e-01
   1.3474e-01
   3.3771e-02
   1.1300e-02
   1.2110e-02
   1.8537e-02
   7.0313e-03
   5.9189e-03
   3.9103e-03
   4.4253e-03
   5.1541e-03
   6.0097e-04
  -2.7835e-05
   1.4314e-03
   1.6354e-03
   8.6584e-04
   4.8698e-04
  -8.6823e-06
   4.7618e-04]'; %19, 7.7391
OptDat.Pin{19}  = [7.3139e+00
   5.3160e+00
   1.0234e+01
   3.0118e-01
   1.3460e-01
   3.3565e-02
   1.1134e-02
   1.2150e-02
   1.7883e-02
   7.0463e-03
   5.9680e-03
   3.8603e-03
   4.4734e-03
   4.9223e-03
   6.0021e-04
  -1.5100e-04
   1.3716e-03
   1.5212e-03
   9.2927e-04
   5.6634e-04
   1.1634e-09]'; %18, 7.7625
OptDat.Pin{18}  = [7.3137e+00
   5.3158e+00
   1.0235e+01
   3.0106e-01
   1.3456e-01
   3.3567e-02
   1.1217e-02
   1.2147e-02
   1.7871e-02
   7.0895e-03
   6.0126e-03
   3.9083e-03
   4.4611e-03
   4.8941e-03
   6.5091e-04
  -1.1394e-04
   1.3705e-03
   1.5247e-03
   9.2963e-04
   6.1535e-04]'; %17, 7.7616
OptDat.Pin{17}  = [7.3157e+00
   5.3175e+00
   1.0231e+01
   3.0268e-01
   1.3512e-01
   3.2237e-02
   1.0909e-02
   1.2331e-02
   1.7526e-02
   6.6793e-03
   5.4193e-03
   3.5147e-03
   4.8426e-03
   4.5611e-03
   4.9343e-04
  -4.2143e-04
   1.0555e-03
   1.6787e-03
   8.5945e-04]'; %16, 7.7832
OptDat.Pin{16}  = [7.3152e+00
   5.3180e+00
   1.0231e+01
   3.0361e-01
   1.3529e-01
   3.1431e-02
   1.0516e-02
   1.1832e-02
   1.7493e-02
   6.7895e-03
   4.8870e-03
   3.0067e-03
   4.4872e-03
   4.5746e-03
   6.2103e-05
  -5.5232e-04
   1.0385e-03
   1.7401e-03]'; %15, 7.8184
OptDat.Pin{15}  = [7.3087e+00
   5.3118e+00
   1.0243e+01
   3.0473e-01
   1.3158e-01
   3.1929e-02
   1.1613e-02
   1.1348e-02
   1.7660e-02
   5.7205e-03
   5.0414e-03
   3.6746e-03
   3.7709e-03
   3.9684e-03
   4.9322e-05
  -3.8995e-04
   2.0908e-03]'; %14, 7.8854
OptDat.Pin{14}  = [7.3131e+00
   5.3155e+00
   1.0236e+01
   3.0462e-01
   1.3026e-01
   3.0789e-02
   1.0527e-02
   1.2112e-02
   1.6399e-02
   5.6039e-03
   4.4397e-03
   3.1699e-03
   3.4349e-03
   3.4021e-03
  -1.3638e-05
  -4.7074e-05]'; %13, 7.9648
OptDat.Pin{13}  = [7.3137e+00
   5.3158e+00
   1.0235e+01
   3.0485e-01
   1.3071e-01
   3.1851e-02
   1.0183e-02
   1.1552e-02
   1.6246e-02
   5.7328e-03
   4.2017e-03
   3.0164e-03
   3.0013e-03
   3.4854e-03
  -1.8741e-04]'; %12, 7.9637
OptDat.Pin{12}  = [7.3169e+00
   5.3177e+00
   1.0232e+01
   3.0434e-01
   1.3098e-01
   3.2026e-02
   9.8591e-03
   1.1903e-02
   1.6018e-02
   5.9079e-03
   4.4849e-03
   3.0438e-03
   3.0857e-03
   3.3798e-03]'; %11, 7.9629
OptDat.Pin{11}  = [7.3069e+00
   5.3087e+00
   1.0249e+01
   3.0907e-01
   1.3196e-01
   2.4946e-02
   1.1516e-02
   1.0794e-02
   1.5381e-02
   4.1411e-03
   4.6029e-03
   3.6014e-03
   3.6383e-03]'; %10, 8.0364
OptDat.Pin{10}  = [7.3036e+00
   5.3068e+00
   1.0253e+01
   3.0735e-01
   1.2922e-01
   3.0190e-02
   1.2564e-02
   8.1778e-03
   1.2092e-02
   2.3942e-03
   7.1654e-03
   4.2732e-03]'; %9, 8.1108
OptDat.Pin{9}  = [7.3188e+00
   5.3209e+00
   1.0225e+01
   3.1040e-01
   1.2845e-01
   2.5425e-02
   1.1363e-02
   8.8466e-03
   9.3977e-03
   2.8727e-03
   6.9662e-03]'; %8, 8.2334
OptDat.Pin{8}  = [7.3392e+00
   5.3159e+00
   1.0238e+01
   3.1288e-01
   1.2470e-01
   2.4597e-02
   7.0623e-03
   9.6340e-03
   1.2487e-02
   3.6397e-03]'; %7, 8.3369
OptDat.Pin{7}  = [7.2973e+00
   5.2983e+00
   1.0267e+01
   3.1589e-01
   1.1868e-01
   1.9869e-02
   9.3025e-03
   8.8679e-03
   1.3574e-02]'; %6, 8.3617
OptDat.Pin{6}  = [7.3738e+00
   5.3618e+00
   1.0134e+01
   3.1352e-01
   1.2566e-01
   1.6936e-02
   3.6205e-03
   1.4348e-02]'; %5, 8.7478
OptDat.Pin{5}  = [7.3145e+00
   5.3167e+00
   1.0233e+01
   3.1729e-01
   1.1305e-01
   1.3167e-02
   9.1453e-03]'; %4, 8.9271
OptDat.Pin{4}  = [7.3264e+00
   5.3264e+00
   1.0213e+01
   3.2016e-01
   1.1159e-01
   1.1105e-02]'; %3, 8.9652
OptDat.Pin{3}  = [7.3250e+00
   5.3253e+00
   1.0216e+01
   3.2605e-01
   1.1057e-01]'; %2, 9.0127
OptDat.Pin{2}  = [7.3250e+00
   5.3253e+00
   1.0216e+01
   3.2605e-01]'; %1, 1.1388e+01
OptDat.Pin{1}  = [7.4213e+00
   5.3916e+00
   1.0057e+01]'; %0, 2.7621e+01


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
