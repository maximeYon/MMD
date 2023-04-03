clear all

spec_fn = '/Users/daniel/Dropbox/NMRdata/AV4_500/Simon/MIC5_kasia/47/Itd1Dat.mat';

load(spec_fn)

maxI = max(Itd1(:));
thresh = 0.0002*maxI;

mask = any(Itd1>thresh,2);
mask_ind = find(mask);

cmin = thresh;
cmax = .02*maxI;

im2d = Itd1';
im2d(im2d<cmin) = cmin;
im2d(im2d>cmax) = cmax;
im2d = log10(im2d);

figure(1), clf
imagesc(im2d)
set(gca,'YDir','normal')
colormap('hot')

opt = mdm_opt();
opt = ilt_opt(opt);
opt.ilt.ind_start = 1;
opt.ilt.n_in = 100;
opt.ilt.n_out = 10;
opt.ilt.n_proliferation = 2;
opt.ilt.n_extinction = 2;
opt.ilt.dmin = 1e-12;
opt.ilt.dmax = .2e-8;

xps.b = b;
xps.n = numel(b);
xps.b_delta = ones(xps.n,1);

Imasked = Itd1(mask,:);

Nchannels = sum(mask,1);
mmasked = zeros(Nchannels,2*opt.ilt.n_out+1);

Nx = 100;
dist_s.x = linspace(log10(opt.ilt.dmin),log10(opt.ilt.dmax),Nx)';
dist_s.xsigma = 3*(dist_s.x(2) - dist_s.x(1));
distmasked = zeros(Nchannels,Nx);

tic
for nchannel = 1:Nchannels
    s = Imasked(nchannel,:)';
    
    ind = 1:xps.n;
    m = ilt_1d_data2fit(s, xps, opt, ind);
    mmasked(nchannel,:) = m;
    
    dd = ilt_m2dd(m);
    [dist_d.n,D,dist_d.w] = ilt_dist2par(dd);
    dist_d.x = log10(D);
    dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
    distmasked(nchannel,:) = dist_s.w;
    
    s_fit = ilt_1d_fit2data(m, xps);
        
%     figure(2), clf, subplot(1,2,1), semilogy(b,s,'o',b,s_fit,'-'), title(num2str(mask_ind(nchannel))), subplot(1,2,2), plot(dist_s.x,dist_s.w,'-'), pause(.1)
end
toc

%%
dist = zeros(size(Itd1,1),Nx);
dist(mask_ind,:) = distmasked;

maxdist = max(dist(:));
thresh = 0.0002*maxdist;

cmin = thresh;
cmax = .02*maxdist;

im2d = dist';
im2d(im2d<cmin) = cmin;
im2d(im2d>cmax) = cmax;
im2d = log10(im2d);

im1d_x = sum(dist',1);
im1d_y = sum(dist',2);

figure(2), clf
axh_2d = axes('position',[.3 .1 .6 .6]);
imagesc(ppm,dist_s.x,im2d)
set(axh_2d,'YDir','normal','XDir','reverse','YAxisLocation','right')
colormap('hot')

axh_projx = axes('position',[.3 .75 .6 .2]);
plot(ppm,im1d_x,'-')
axis('tight')
set(axh_projx,'XDir','reverse','YLim',.05*max(im1d_x)*[-.1 1.1])

axh_projy = axes('position',[.02 .1 .2 .6]);
plot(im1d_y,dist_s.x,'-')
axis('tight')
set(axh_projy,'XDir','reverse','XLim',.1*max(im1d_y)*[-.1 1.1],'YLim',[min(dist_s.x) max(dist_s.x)])


