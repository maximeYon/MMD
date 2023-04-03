clear all

wd = cd;
fpath = wd;

DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
ExpNam = 'C14E5_5'; expno = 16;

fpath = [DataDir '/' ExpNam '/' num2str(expno)];

Gmax = 0.3; tau = 40e-3;
epsilon_up = .1; plateau = .0001; epsilon_down = .15;

Gmax = 0.3; tau = 25e-3;
epsilon_up = .1; plateau = .0001; epsilon_down = .15;

Np = 1000; %number of points in waveform

%timing parameters relative to a total echo time = 1
Deltapsi = 2*pi; zeta = acos(1/sqrt(3));
%Deltapsi = 1*pi; zeta = pi/2;
%Deltapsi = 2*pi; zeta = pi/2;
theta = 0; phi = 0;

%------------------------
tE = 1;
dt = tE/Np;
t = tE*linspace(0,1,Np);

Np_epsilon_up = round(epsilon_up/dt);
t_epsilon_up = pi/2*linspace(0,1,Np_epsilon_up)';
G_up = sin(t_epsilon_up);
%figure(1), clf, plot(t_epsilon_up,G_up,'-'), return

Np_epsilon_down = round(epsilon_down/dt);
t_epsilon_down = pi/2*linspace(-1,1,Np_epsilon_down)';
G_down = 1 - .5*(1+sin(t_epsilon_down));
%figure(1), clf, plot(t_epsilon_down,G_down,'-'), return

Np_plateau = round(plateau/dt);
Np_inter = round(tE/dt)-2*(Np_epsilon_up+Np_plateau+Np_epsilon_down);

G = [G_up; ones(Np_plateau,1); G_down; zeros(Np_inter,1); G_down-1; -1*ones(Np_plateau,1); -flipud(G_up)];

%figure(1), clf, plot(t,G,'-'), return


psi0 = 0;

F = cumsum(G);
F = F/max(F);

td = sum(F.^2*dt);
psirate = Deltapsi/td*F.^2;
psi = psi0 + cumsum(psirate*dt);

qx = F.*sin(zeta).*cos(psi);
qy = F.*sin(zeta).*sin(psi);
qz = F.*cos(zeta);

qx_old = qx;
qz_old = qz;
qx = qx_old*cos(theta) + qz_old*sin(theta);
qz = qz_old*cos(theta) - qx_old*sin(theta);
qx_old = qx;
qy_old = qy;
qx = qx_old*cos(phi) - qy_old*sin(phi);
qy = qy_old*cos(phi) + qx_old*sin(phi);


Giso_x = gradient(qx)/dt;
Giso_y = gradient(qy)/dt;
Giso_z = gradient(qz)/dt;
Gdir = gradient(F)/dt;

Giso_x = [0; diff(qx)/dt];
Giso_y = [0; diff(qy)/dt];
Giso_z = [0; diff(qz)/dt];
Gdir = [0; diff(F)/dt];

Giso_m = sqrt(Giso_x.^2 + Giso_y.^2 + Giso_z.^2);

maxG = max(abs([Giso_x; Giso_y; Giso_z; Giso_m; Gdir]));

Giso_x = Giso_x/maxG;
Giso_y = Giso_y/maxG;
Giso_z = Giso_z/maxG;
Giso_m = Giso_m/maxG;
Gdir = Gdir/maxG;

qx = cumsum(Giso_x*dt);
qy = cumsum(Giso_y*dt);
qz = cumsum(Giso_z*dt);
qm = sqrt(qx.^2 + qy.^2 + qz.^2);
qdir = cumsum(Gdir*dt);

b = sum(qm.^2*dt);

bxx = sum(qx.^2*dt);
byy = sum(qy.^2*dt);
bzz = sum(qz.^2*dt);
bxy = sum(qx.*qy*dt);
bxz = sum(qx.*qz*dt);
byz = sum(qy.*qz*dt);

gamma = 26.75e7;
dt = tau/Np;
t = tau*linspace(0,1,Np)';
qdir = gamma*Gmax*cumsum(Gdir*dt);
b = sum(qdir.^2*dt)
C = b/(gamma^2*Gmax^2*tE^3)

dGxdt = gradient(Gmax*Giso_x/(dt));
dGydt = (Gmax*gradient(Giso_y)/dt);
dGzdt = (Gmax*gradient(Giso_z)/dt);
dGmdt = (Gmax*gradient(Giso_m)/dt);
dGdirdt = (Gmax*gradient(Gdir)/dt);

figure(1), clf
subplot(2,1,1)
plot(t,Gmax*Giso_x,'r-',t,Gmax*Giso_y,'g-',t,Gmax*Giso_z,'b-',t,Gmax*Giso_m,'k-.',t,Gmax*Gdir,'k-')
ylabel('G / Tm^-^1')
title(['b = ' num2str(b/1e9,2) ' 10^9 sm^-^2'])

subplot(2,1,2)
plot(t,dGxdt,'r-',t,dGydt,'g-',t,dGzdt,'b-',t,dGmdt,'k-.',t,dGdirdt,'k-')
xlabel('t / s'), ylabel('(dG/dt) / Tm^-^1s^-^1')

td = sum(F.^2*dt)
q = max(qdir)
b=q^2*td

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


fprintf(fidx,'%s\n',['##NPOINTS=' num2str(Np)]);
fprintf(fidx,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidy,'%s\n',['##NPOINTS=' num2str(Np)]);
fprintf(fidy,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidz,'%s\n',['##NPOINTS=' num2str(Np)]);
fprintf(fidz,'%s\n',['##XYDATA= (X++(Y..Y))']);
fprintf(fidr,'%s\n',['##NPOINTS=' num2str(Np)]);
fprintf(fidr,'%s\n',['##XYDATA= (X++(Y..Y))']);

for nG = 1:length(Giso_x)
    fprintf(fidx,'%f\n',Giso_x(nG));
    fprintf(fidxtxt,'%f\n',Giso_x(nG));
    fprintf(fidy,'%f\n',Giso_y(nG));
    fprintf(fidytxt,'%f\n',Giso_y(nG));
    fprintf(fidz,'%f\n',Giso_z(nG));
    fprintf(fidztxt,'%f\n',Giso_z(nG));
    fprintf(fidr,'%f\n',Gdir(nG));
    fprintf(fidrtxt,'%f\n',Gdir(nG));
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
