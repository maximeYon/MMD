function bruker_axde_protocol(run_fn)
% function bruker_axde_protocol(run_fn)
%
% run_fn - where to put the data
%
% Acquisition protocol for axisymmetric diffusion encoding.
%
% Used in Topgaard, Phys. Chem. Chem. Phys. 18, 8545 (2016).
% http://dx.doi.org/10.1039/c5cp07251d

clear all;

% Define path for output folders
[run_path,run_name,run_ext] = fileparts(run_fn);
out_path = fullfile(run_path,'protocol');
framework_path = fileparts(fileparts(run_path));
pa_path = fullfile(framework_path,'tools','uvec','repulsion_angles');


% Define timing parameters relative to a total echo time = 1
n_b = 4;
bmin = .01;
%brel = linspace(bmin,1,n_b)'; %linear spacing b
brel = logspace(log10(bmin),0,n_b)'; %linear spacing b
n_bd = 4; %n_bd = 4, 7, 10, 13,...
xps.b_delta = linspace(1,-.5,n_bd)';
%xps.b_delta = [1 0 -.5 -sqrt(.125) sqrt(.125) .5 sqrt(.5) sqrt(.75)]'; n_bd = numel(xps.b_delta);
n_br = 11;
load(fullfile(pa_path,num2str(n_br)))
%xps.b_delta = [0]; n_bd = numel(xps.b_delta);
%theta = [0]; phi = [0]; n_br = numel(theta);

[b_ind,bd_ind,br_ind] = ndgrid(1:n_b,1:n_bd,1:n_br);
xps.n = numel(b_ind);

xps.b_ind = reshape(b_ind,xps.n,1); % Same b_ind means same b-tensor size
xps.bd_ind = reshape(bd_ind,xps.n,1); % Same bd_ind means same b-tensor shape
xps.br_ind = reshape(br_ind,xps.n,1); % Same br_ind means same b-tensor orientation
[a_ind,~] = ndgrid(1:(n_b*n_bd),1:n_br);
xps.a_ind = reshape(a_ind,xps.n,1); % Same a_ind means same b-tensor size and shape

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
subplot(2,1,1)
hx = plot(1:xps.n,g.xa,'ro',1:xps.n,g.xb,'rs',1:xps.n,g.xc,'rx');
hold on
hy = plot(1:xps.n,g.ya,'go',1:xps.n,g.yb,'gs',1:xps.n,g.yc,'gx');
hz = plot(1:xps.n,g.za,'bo',1:xps.n,g.zb,'bs',1:xps.n,g.zc,'bx');
set(hx,'MarkerSize',10), set(hy,'MarkerSize',8), set(hz,'MarkerSize',6)
title(['td1 = ' num2str(xps.n)])
subplot(2,1,2)
ha = plot(1:xps.n,sqrt(g.xa.^2+g.ya.^2+g.za.^2),'ro');
hb = plot(1:xps.n,sqrt(g.xa.^2+g.ya.^2+g.za.^2),'ro');
hc = plot(1:xps.n,sqrt(g.xa.^2+g.ya.^2+g.za.^2),'ro');

param = {'xa','xb','xc','ya','yb','yc','za','zb','zc'};
for nparam = 1:numel(param)
    eval(['gtemp = g.' param{nparam} ';'])
    out_fn = fullfile(out_path,['g' param{nparam}]);
    res = bruker_mkshapefile(out_fn,gtemp,param{nparam});
end

save(fullfile(out_path,'g_xps'), 'xps')
