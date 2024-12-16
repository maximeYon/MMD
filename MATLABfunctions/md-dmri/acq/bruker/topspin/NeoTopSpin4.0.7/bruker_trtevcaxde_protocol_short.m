% Acquisition protocol for axisymmetric diffusion encoding.
%
% Used in Topgaard, Phys. Chem. Chem. Phys. 18, 8545 (2016).
% http://dx.doi.org/10.1039/c5cp07251d
%
% Transfer all files from the protocol folder to
% /opt/topspin<version>/exp/stan/nmr/lists/gp/user
% and <dataset>/<expno>

clear all;

% Define path for output folders
%[run_path,run_name,run_ext] = fileparts(run_fn);
run_path = fileparts(mfilename('fullpath'));
%out_path = fullfile(run_path,'protocol_t2axde');
out_path = '/opt/topspin4.0.7/exp/stan/nmr/lists/gp/user';
%out_path = '/opt/topspin3.5pl4/exp/stan/nmr/lists/gp/user';
%out_path = '/opt/topspin3.2/exp/stan/nmr/lists/gp/user';

framework_path = fileparts(fileparts(fileparts(fileparts(run_path))));
pa_path = fullfile(framework_path,'tools','uvec','repulsion_angles');

% Define timing parameters relative to a total echo time = 1
n_b = 2;
bmin = 0.001;
%brel = linspace(bmin,1,n_b)'; %linear spacing b
%brel = (linspace(sqrt(bmin),1,n_b)').^2; %quadratic spacing b
brel = (linspace(sqrt(bmin),1,n_b+1)').^2; brel = brel(2:(n_b+1));%quadratic spacing
% brel = logspace(log10(bmin),0,n_b)'; %logarithmic spacing b
n_bd = 4; %n_bd = 4, 7, 10, 13,...
%xps.b_delta = linspace(-.5,1,n_bd)';
%xps.b_delta = [1 0 -.5 -sqrt(.125) sqrt(.125) .5 sqrt(.5) sqrt(.75)]'; n_bd = numel(xps.b_delta);
%xps.b_delta = [1 0 -.5 .5 sqrt(.5) sqrt(.75)]'; n_bd = numel(xps.b_delta);
xps.b_delta = [-.5 .5 1]'; n_bd = numel(xps.b_delta);
%xps.b_delta = [-.5 0 1]'; n_bd = numel(xps.b_delta);
%xps.b_delta = [0]; n_bd = numel(xps.b_delta);
n_br = 23;
load(fullfile(pa_path,num2str(n_br)))
%theta = pi/2*[0]; phi = pi/2*[0]; n_br = numel(theta);
%theta = pi/2*[0 2 1 1 1 1 .5 1.5 1 1 .5 1.5]; phi = pi/2*[0 0 0 2 1 3 0 0 .5 1.5 1 3]; n_br = numel(theta);
%theta = 2*pi*linspace(0,1,17); phi = 3*pi/4*ones(size(theta)); n_br = numel(theta);
%phi = 2*pi*linspace(0,1,9); theta = 1*pi/4*ones(size(phi)); n_br = numel(theta);
%theta = 2*pi*linspace(0,1,33); phi = 2*theta; n_br = numel(theta);
%theta = acos(2*rand(n_br,1)-1); phi = 2*pi*rand(n_br,1);

%[b_ind,bd_ind,br_ind,te_ind,tr_ind] = ndgrid(1:n_b,1:n_bd,1:n_br,1:n_te,1:n_tr);
%[b_ind,bd_ind,br_ind] = ndgrid(1:n_b,1:n_bd,1:n_br);
[b_ind,bd_ind,br_ind] = ndgrid(1:n_b,1:n_bd,1:n_br);
xps.n = numel(b_ind);

xps.b_ind = reshape(b_ind,xps.n,1); % Same b_ind means same b-tensor size
xps.bd_ind = reshape(bd_ind,xps.n,1); % Same bd_ind means same b-tensor shape
xps.br_ind = reshape(br_ind,xps.n,1); % Same br_ind means same b-tensor orientation
% xps.te_ind = reshape(te_ind,xps.n,1); % Same te_ind means same echo time
% xps.tr_ind = reshape(tr_ind,xps.n,1); % Same tr_ind means same recovery time
% [a_ind,~] = ndgrid(1:(n_b*n_bd),1:n_br,1:n_te,1:n_tr);
[a_ind,~] = ndgrid(1:(n_b*n_bd),1:n_br);
xps.a_ind = reshape(a_ind,xps.n,1); % Same a_ind means same b-tensor size and shape

brel = brel(xps.b_ind);
xps.b_delta = xps.b_delta(xps.bd_ind);
alpha = phi(xps.br_ind);
beta = theta(xps.br_ind);
gamma = zeros(size(xps.br_ind));
%figure(1), clf, plot(alpha,beta,'o'), return

%n_bs = 16*n_b-1;
n_bs = 32;
brel = cat(1,(linspace(sqrt(bmin),1,n_bs)').^2,brel);
xps.b_delta = cat(1,zeros(n_bs,1),xps.b_delta);
alpha = cat(1,zeros(n_bs,1),alpha);
beta = cat(1,zeros(n_bs,1),beta);
gamma = cat(1,zeros(n_bs,1),gamma);
xps.n = numel(brel);

grel = sqrt(brel);
zeta = acos(sqrt((xps.b_delta*2+1)/3));
g.a = sin(zeta).*grel;
g.b = sin(zeta).*grel;
g.c = cos(zeta).*grel;


param = {'xx','xy','xz','yx','yy','yz','zx','zy','zz'};
for nparam = 1:numel(param)
    eval(['R.' param{nparam} ' = zeros(xps.n,1);'])
end

for n = 1:xps.n
    [rotmat,~] = tm_euler_angles2rotmat(alpha(n),beta(n),gamma(n));
    
    R.xx(n,1) = rotmat(1,1);
    R.xy(n,1) = rotmat(1,2);
    R.xz(n,1) = rotmat(1,3);
    R.yx(n,1) = rotmat(2,1);
    R.yy(n,1) = rotmat(2,2);
    R.yz(n,1) = rotmat(2,3);
    R.zx(n,1) = rotmat(3,1);
    R.zy(n,1) = rotmat(3,2);
    R.zz(n,1) = rotmat(3,3);
end

nblock = 1;
n_tr = 32;
vdte = 1e-3;
vdtr = flipud(logspace(log10(10e-3),log10(15),n_tr)');
%vdtr = flipud(logspace(log10(10e-3),log10(1),n_tr)');
vc = 0;

    sv_cell{nblock}.g.xa = zeros(numel(vdtr),1);
    sv_cell{nblock}.g.xb = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.xc = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.ya = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.yb = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.yc = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.za = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.zb = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.zc = sv_cell{nblock}.g.xa;

    sv_cell{nblock}.vdte = vdte*ones(numel(sv_cell{nblock}.g.xa),1);
    sv_cell{nblock}.vdtr = vdtr;
    sv_cell{nblock}.vc = vc*ones(numel(sv_cell{nblock}.g.xa),1);

nblock = 2;
n_te = 32;
vdte = logspace(log10(1e-3),log10(200e-3),n_te)';
vdtr = 10;
vc = 0;

    sv_cell{nblock}.g.xa = zeros(numel(vdte),1);
    sv_cell{nblock}.g.xb = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.xc = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.ya = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.yb = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.yc = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.za = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.zb = sv_cell{nblock}.g.xa;
    sv_cell{nblock}.g.zc = sv_cell{nblock}.g.xa;

    sv_cell{nblock}.vdte = vdte;
    sv_cell{nblock}.vdtr = vdtr*ones(numel(sv_cell{nblock}.g.xa),1);
    sv_cell{nblock}.vc = vc*ones(numel(sv_cell{nblock}.g.xa),1);

nblock = 3;
n_tr = 32;
vdte = 1e-3;
vdtr = flipud(logspace(log10(10e-3),log10(15),n_tr)');
vc = 1;
brel = .25;
b_delta = 0;

grel = sqrt(brel)*ones(size(vdtr));
zeta = acos(sqrt((b_delta*2+1)/3));

    sv_cell{nblock}.g.xa = sin(zeta).*grel;
    sv_cell{nblock}.g.xb = zeros(size(grel));
    sv_cell{nblock}.g.xc = zeros(size(grel));
    sv_cell{nblock}.g.ya = zeros(size(grel));
    sv_cell{nblock}.g.yb = sin(zeta).*grel;
    sv_cell{nblock}.g.yc = zeros(size(grel));
    sv_cell{nblock}.g.za = zeros(size(grel));
    sv_cell{nblock}.g.zb = zeros(size(grel));
    sv_cell{nblock}.g.zc = cos(zeta).*grel;
    
sv_cell{nblock}.vdte = vdte*ones(numel(sv_cell{nblock}.g.xa),1);
sv_cell{nblock}.vdtr = vdtr;
sv_cell{nblock}.vc = vc*ones(numel(sv_cell{nblock}.g.xa),1);

    
 
nblock = 4;
n_te = 32;
vdte = logspace(log10(1e-3),log10(200e-3),n_te)';
vdtr = 10;
vc = 1;
brel = .25;
b_delta = 0;

grel = sqrt(brel)*ones(size(vdte));
zeta = acos(sqrt((b_delta*2+1)/3));

    sv_cell{nblock}.g.xa = sin(zeta).*grel;
    sv_cell{nblock}.g.xb = zeros(size(grel));
    sv_cell{nblock}.g.xc = zeros(size(grel));
    sv_cell{nblock}.g.ya = zeros(size(grel));
    sv_cell{nblock}.g.yb = sin(zeta).*grel;
    sv_cell{nblock}.g.yc = zeros(size(grel));
    sv_cell{nblock}.g.za = zeros(size(grel));
    sv_cell{nblock}.g.zb = zeros(size(grel));
    sv_cell{nblock}.g.zc = cos(zeta).*grel;

    sv_cell{nblock}.vdte = vdte;
    sv_cell{nblock}.vdtr = vdtr*ones(numel(sv_cell{nblock}.g.xa),1);
    sv_cell{nblock}.vc = vc*ones(numel(sv_cell{nblock}.g.xa),1);
    
% nblock = 3;
% vdte = 1e-3;
% vdtr = 5;
% vc = 1;
% bmin = 0.01;
% n_b = 25;
% brel = (linspace(sqrt(bmin),1,n_b)').^2;
% %brel = logspace(log10(bmin),0,n_b)'; %logarithmic spacing b
% b_delta = 0;
% grel = sqrt(brel);
% zeta = acos(sqrt((b_delta*2+1)/3));
% 
%     sv_cell{nblock}.g.xa = sin(zeta).*grel;
%     sv_cell{nblock}.g.xb = zeros(n_b,1);
%     sv_cell{nblock}.g.xc = zeros(n_b,1);
%     sv_cell{nblock}.g.ya = zeros(n_b,1);
%     sv_cell{nblock}.g.yb = sin(zeta).*grel;
%     sv_cell{nblock}.g.yc = zeros(n_b,1);
%     sv_cell{nblock}.g.za = zeros(n_b,1);
%     sv_cell{nblock}.g.zb = zeros(n_b,1);
%     sv_cell{nblock}.g.zc = cos(zeta).*grel;
% 
%     sv_cell{nblock}.vdte = vdte*ones(numel(sv_cell{nblock}.g.xa),1);
%     sv_cell{nblock}.vdtr = vdtr*ones(numel(sv_cell{nblock}.g.xa),1);
%     sv_cell{nblock}.vc = vc*ones(numel(sv_cell{nblock}.g.xa),1);

    
% vdte = [1          1   1 30 100         30 ]*1e-3;
% vdtr = 5*[1 .1*sqrt(10) .1  1   1 .1*sqrt(10)];
vdte = [1 25  100  1  1  25]*1e-3;
vdtr = [10 10   10   2 .5 2];
% vdte = [1]*1e-3;
% vdtr = 5*[1];
vc = ones(size(vdte));

param = {'xa','xb','xc','ya','yb','yc','za','zb','zc'};
param_v = {'vdtr','vdte','vc'};

for nvd = 1:numel(vdte)
    nblock = nblock+1;
    sv_cell{nblock}.g.xa = R.xx.*g.a;
    sv_cell{nblock}.g.xb = R.xy.*g.b;
    sv_cell{nblock}.g.xc = R.xz.*g.c;
    sv_cell{nblock}.g.ya = R.yx.*g.a;
    sv_cell{nblock}.g.yb = R.yy.*g.b;
    sv_cell{nblock}.g.yc = R.yz.*g.c;
    sv_cell{nblock}.g.za = R.zx.*g.a;
    sv_cell{nblock}.g.zb = R.zy.*g.b;
    sv_cell{nblock}.g.zc = R.zz.*g.c;

    sv_cell{nblock}.vdte = vdte(nvd)*ones(numel(sv_cell{nblock}.g.xa),1);
    sv_cell{nblock}.vdtr = vdtr(nvd)*ones(numel(sv_cell{nblock}.g.xa),1);
    sv_cell{nblock}.vc = vc(nvd)*ones(numel(sv_cell{nblock}.g.xa),1);
    
    for nblock1 = 1:4
        nblock = nblock+1;
        for nparam = 1:numel(param)
            sv_cell{nblock}.g.(param{nparam}) = sv_cell{nblock1}.g.(param{nparam});
        end
        for nparam = 1:numel(param_v)
            sv_cell{nblock}.(param_v{nparam}) = sv_cell{nblock1}.(param_v{nparam});
        end
    end
   
end

clear g
for nparam = 1:numel(param)
    g.(param{nparam}) = [];
end
for nparam = 1:numel(param_v)
    sv.(param_v{nparam}) = [];
end


for nblock = 1:numel(sv_cell)
    for nparam = 1:numel(param)
        g.(param{nparam}) = cat(1,g.(param{nparam}),sv_cell{nblock}.g.(param{nparam}));
    end
    for nparam = 1:numel(param_v)
        sv.(param_v{nparam}) = cat(1,sv.(param_v{nparam}),sv_cell{nblock}.(param_v{nparam}));
    end
end

xps.n = numel(g.xa);

figure(11), clf
subplot(3,1,1)
hx = plot(1:xps.n,g.xa,'ro',1:xps.n,g.xb,'rs',1:xps.n,g.xc,'rx');
hold on
hy = plot(1:xps.n,g.ya,'go',1:xps.n,g.yb,'gs',1:xps.n,g.yc,'gx');
hz = plot(1:xps.n,g.za,'bo',1:xps.n,g.zb,'bs',1:xps.n,g.zc,'bx');
set(hx,'MarkerSize',10), set(hy,'MarkerSize',8), set(hz,'MarkerSize',6)
title(['td1 = ' num2str(xps.n)])
subplot(3,1,2)
ha = plot(1:xps.n,sqrt(g.xa.^2+g.ya.^2+g.za.^2),'ro');
hold on
hb = plot(1:xps.n,sqrt(g.xb.^2+g.yb.^2+g.zb.^2),'gs');
hc = plot(1:xps.n,sqrt(g.xc.^2+g.yc.^2+g.zc.^2),'bx');
subplot(3,1,3)
hte = semilogy(1:xps.n,sv.vdte,'ro',1:xps.n,sv.vdtr,'bo');

param = {'xa','xb','xc','ya','yb','yc','za','zb','zc'};
for nparam = 1:numel(param)
    eval(['gtemp = g.' param{nparam} ';'])
    out_fn = fullfile(out_path,['g' param{nparam}]);
    res = bruker_mkshapefile(out_fn,gtemp,param{nparam});
end


td1 = xps.n;

lists_path = '/opt/topspin4.0.7/exp/stan/nmr/lists/';

fid_tr = fopen(fullfile(lists_path,'vd','DT_vdtr'),'w');
fid_te = fopen(fullfile(lists_path,'vd','DT_vdte'),'w');
fid_vc = fopen(fullfile(lists_path,'vc','DT_vc'),'w');

for n = 1:td1
    fprintf(fid_tr,'%1.8f\n',sv.vdtr(n));
    fprintf(fid_te,'%1.8f\n',sv.vdte(n));
    fprintf(fid_vc,'%i\n',sv.vc(n));
end

fclose(fid_tr);
fclose(fid_te);
fclose(fid_vc);