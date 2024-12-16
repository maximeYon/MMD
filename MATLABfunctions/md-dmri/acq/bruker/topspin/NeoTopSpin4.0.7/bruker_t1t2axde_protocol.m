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
n_b = 8;
bmin = 0.01;
%brel = linspace(bmin,1,n_b)'; %linear spacing b
brel = (linspace(sqrt(bmin),1,n_b)').^2; %quadratic spacing b
% brel = logspace(log10(bmin),0,n_b)'; %logarithmic spacing b
n_bd = 4; %n_bd = 4, 7, 10, 13,...
xps.b_delta = linspace(1,-.5,n_bd)';
%xps.b_delta = [1 0 -.5 -sqrt(.125) sqrt(.125) .5 sqrt(.5) sqrt(.75)]'; n_bd = numel(xps.b_delta);
%xps.b_delta = [1 0 -.5 .5 sqrt(.5) sqrt(.75)]'; n_bd = numel(xps.b_delta);
%xps.b_delta = [0 1]'; n_bd = numel(xps.b_delta);
%xps.b_delta = [0]; n_bd = numel(xps.b_delta);
n_br = 15;
load(fullfile(pa_path,num2str(n_br)))
%theta = pi/2*[0]; phi = pi/2*[0]; n_br = numel(theta);
%theta = pi/2*[0 2 1 1 1 1 .5 1.5 1 1 .5 1.5]; phi = pi/2*[0 0 0 2 1 3 0 0 .5 1.5 1 3]; n_br = numel(theta);
%theta = 2*pi*linspace(0,1,17); phi = 3*pi/4*ones(size(theta)); n_br = numel(theta);
%phi = 2*pi*linspace(0,1,9); theta = 1*pi/4*ones(size(phi)); n_br = numel(theta);
%theta = 2*pi*linspace(0,1,33); phi = 2*theta; n_br = numel(theta);
%theta = acos(2*rand(n_br,1)-1); phi = 2*pi*rand(n_br,1);
n_te = 8; %
vdte = logspace(log10(1e-3),log10(200e-3),n_te)';
n_tr = 8; %
vdtr = flipud(logspace(log10(200e-3),log10(5000e-3),n_tr)');

[b_ind,bd_ind,br_ind,te_ind,tr_ind] = ndgrid(1:n_b,1:n_bd,1:n_br,1:n_te,1:n_tr);
xps.n = numel(b_ind);

xps.b_ind = reshape(b_ind,xps.n,1); % Same b_ind means same b-tensor size
xps.bd_ind = reshape(bd_ind,xps.n,1); % Same bd_ind means same b-tensor shape
xps.br_ind = reshape(br_ind,xps.n,1); % Same br_ind means same b-tensor orientation
xps.te_ind = reshape(te_ind,xps.n,1); % Same te_ind means same echo time
xps.tr_ind = reshape(tr_ind,xps.n,1); % Same tr_ind means same recovery time
[a_ind,~] = ndgrid(1:(n_b*n_bd),1:n_br,1:n_te,1:n_tr);
xps.a_ind = reshape(a_ind,xps.n,1); % Same a_ind means same b-tensor size and shape

sv.vdte = vdte(xps.te_ind);
sv.vdtr = vdtr(xps.tr_ind);
brel = brel(xps.b_ind);
xps.b_delta = xps.b_delta(xps.bd_ind);

grel = sqrt(brel);
zeta = acos(sqrt((xps.b_delta*2+1)/3));
g.a = sin(zeta).*grel;
g.b = sin(zeta).*grel;
g.c = cos(zeta).*grel;

alpha = phi(xps.br_ind);
beta = theta(xps.br_ind);
gamma = zeros(size(xps.br_ind));
%figure(1), clf, plot(alpha,beta,'o'), return

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

g.xa = R.xx.*g.a;
g.xb = R.xy.*g.b;
g.xc = R.xz.*g.c;
g.ya = R.yx.*g.a;
g.yb = R.yy.*g.b;
g.yc = R.yz.*g.c;
g.za = R.zx.*g.a;
g.zb = R.zy.*g.b;
g.zc = R.zz.*g.c;

figure(1), clf
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

ListDir = '/opt/topspin4.0.7/exp/stan/nmr/lists/vd';

fpath = [ListDir];

fname = 'vdtr';
fid1 = fopen([fpath '/' fname],'w');
fname = 'vdte';
fid2 = fopen([fpath '/' fname],'w');

for n = 1:td1
    fprintf(fid1,'%1.8f\n',sv.vdtr(n));
    fprintf(fid2,'%1.8f\n',sv.vdte(n));
end

fclose(fid1);
fclose(fid2);
