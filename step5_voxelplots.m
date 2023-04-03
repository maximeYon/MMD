clearvars;close all; clc;

home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

% Prepare paths
% Define full paths to the nifti files to be analyzed 
% data_path = '20211019_110649_test_sample_Hong_test_sample_Hong_1_1\127';
data_path = '20220204_tumor3\60';

data_path = split(data_path,'\'); data_path = join(data_path,filesep,1); data_path = data_path{1};
nii_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path filesep 'pdata_mdd' filesep 'nii_xps'];
bs_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path filesep 'pdata_mdd' filesep 'dtor1r2d' filesep 'bootstraps'];


% method = 'dtod';omega_v = [180 80]*2*pi; omega_v = linspace(53,150,5)*2*pi;
% method = 'dtor1r2d'; omega_v = [150 50]*2*pi; omega_v = linspace(50,150,5)*2*pi;
method = 'dtr1r2d'; 

pdata_path = fileparts(fileparts(bs_path));
voxel_path = fullfile(pdata_path,'voxel');
pmaps_path = fullfile(pdata_path,'pmaps');

% Connect to data
s.nii_fn = fullfile(nii_path, 'data.nii.gz');
s.mask_fn = fullfile(nii_path, 'data_mask.nii.gz');
s.xps = mdm_xps_load(fullfile(nii_path, 'data_xps.mat'));

opt = mdm_opt();

bsno = msf_getdirno(bs_path); %bsno = bsno(1:2);
[dpar_4d,dperp_4d,theta_4d,phi_4d,d0_4d,rpar_4d,rperp_4d,w_4d] = dtod_bsmerge_initializepars;
m_bsmerge = [];
for nbs = 1:numel(bsno)
    mfs_fn   = fullfile(bs_path,num2str(bsno(nbs)),'mfs.mat');
    if exist(mfs_fn,'file')==2
        mfs = mdm_mfs_load(mfs_fn);
        m = mfs.m;
        m_bsmerge = cat(4,m_bsmerge,m);
        [dpar_temp,dperp_temp,theta_temp,phi_temp,d0_temp,rpar_temp,rperp_temp,w_temp] = dtod_4d_m2pars(m);
        [dpar_4d,dperp_4d,theta_4d,phi_4d,d0_4d,rpar_4d,rperp_4d,w_4d] = dtod_bsmerge_catpars(dpar_4d,dperp_4d,theta_4d,phi_4d,d0_4d,rpar_4d,rperp_4d,w_4d,dpar_temp,dperp_temp,theta_temp,phi_temp,d0_temp,rpar_temp,rperp_temp,w_temp);
    end
end


%%

map_nams = {'s0'; 'fractions'; 'mdiso'; 'msddelta'; 'dmdisodnu';'mdii_bin1'; 'dmdisodnu_bin2'};
map_nams = {'s0'};
% map_nams = {'dmdisodnu_bin2'};
Nmaps = numel(map_nams);
for nmap = 1:Nmaps
    map_nam = map_nams{nmap};
    pmap_fn = fullfile(pdata_path,'pmaps',[method '_' map_nam '.nii.gz']);
    [I,nii_h]  = mdm_nii_read(pmap_fn);
    pmaps.(map_nam) = double(I);
end


col_a = [1 0 0
0 1 0
0 0 1
1 .7 0];

nk = 1; % Slice number

% figure(1)
% imagesc(pmaps.s0')
% set(gca,'XDir','reverse','YDir','normal')

% ni_v = [30; 55; 25];% points coordinates
% nj_v = [23; 10; 10];% points coordinates

% ni_v = [30; 52; 23];% points coordinates
% nj_v = [47; 21; 24];% points coordinates

ni_v = [40; 30; 20];% points coordinates
nj_v = [40; 30; 20];% points coordinates

width = 1/Nmaps;

lw = 1;
ms = 2;
fs = 6;
lw_contour = 1;
lw_axes = 1;
lw_plot = .5;
fs_axes = 12;
papersize_i = [5 5];
papersize_p = [10 5*2];
papersize_s = [4 4];
papersize_d = [4 4];
% pformat = '-r600 -dpng';
pformat = '-dpdf';

figure(5), clf
for nmap = 1:Nmaps
    map_nam = map_nams{nmap};
    
    axh = axes('position',[(nmap-1)*width 0 width 1]);

    pmap3d = pmaps.(map_nam);
    if ndims(pmap3d) == 3
        im2d = squeeze(pmap3d(:,:,nk))';
        im2d(sub2ind(size(im2d),nj_v,ni_v))
    elseif ndims(pmap3d) == 4
        im2d = zeros(size(pmap3d,3),size(pmap3d,2),3);
        im2d(:,:,1) = squeeze(pmap3d(1,:,:,nk))'/255;
        im2d(:,:,2) = squeeze(pmap3d(2,:,:,nk))'/255;
        im2d(:,:,3) = squeeze(pmap3d(3,:,:,nk))'/255;
        [im2d(sub2ind(size(im2d),nj_v,ni_v,1*ones(size(ni_v))))...
            im2d(sub2ind(size(im2d),nj_v,ni_v,2*ones(size(ni_v))))...
            im2d(sub2ind(size(im2d),nj_v,ni_v,3*ones(size(ni_v))))]
    elseif ndims(pmap3d) == 2 % mod max
        im2d = pmap3d';
    end
    set(axh,'XLim',[1 size(im2d,1)]);
    set(axh,'YLim',[1 size(im2d,2)]);
    imagesc(axh,im2d)
    oricode = mdm_nii_oricode(nii_h);
    if strcmp(oricode,'LPS')
        set(axh,'YDir','reverse')
    elseif strcmp(oricode,'RAS')
        set(axh,'XDir','reverse','YDir','normal')
    else
        set(axh,'YDir','normal')
    end
 axis(axh,'square','tight','off')
    hold on
    for npixel = 1:numel(ni_v)
        hp = plot(ni_v(npixel),nj_v(npixel),'x','LineWidth',1,'MarkerSize',3);
        set(hp,'Color',col_a(npixel,:))
    end
    
    colormap(axh,'gray')
end

msf_mkdir(voxel_path);

set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize_i],'PaperSize', papersize_i); 
eval(['print ' voxel_path filesep 'labeledvoxels_slice_' num2str(nk) ' -loose ' pformat])

%%
[I,h]  = mdm_nii_read(s.nii_fn);
xps = s.xps;
%%

%Plot distributions for selected voxels

xps_array = [round(xps.b,-7), round(xps.momega,0), round(xps.b_delta,1), xps.theta/pi*180, xps.phi/pi*180];
[xps_array_sort, ind_sort] = sortrows(xps_array,[1 2 3 4 5]);
% ind_sort = 1:xps.n;

figure(5), clf

left = 0.12;
width = .85;
bottom = .08;
dbottom = (1-bottom)/5;
height = dbottom - .02;

axh_b = axes('position',[left bottom+4*dbottom width height]);
ph_b = plot(1:xps.n,xps.b(ind_sort)/1e9,'.');
axh_omega = axes('position',[left bottom+3*dbottom width height]);
ph_omega = plot(1:xps.n,xps.momega(ind_sort)/2/pi,'.');
axh_bd = axes('position',[left bottom+2*dbottom width height]);
ph_bd = plot(1:xps.n,xps.b_delta(ind_sort),'.');
axh_btheta = axes('position',[left bottom+1*dbottom width height]);
ph_btheta = plot(1:xps.n,acos(xps.u((ind_sort),3))/pi*180,'.');
axh_bphi = axes('position',[left bottom+0*dbottom width height]);
ph_bphi = plot(1:xps.n,atan2(xps.u((ind_sort),2),xps.u((ind_sort),1))/pi*180,'.');

set([axh_b; axh_bd; axh_btheta; axh_bphi; axh_omega],...
'XLim',xps.n*[-.05 1.05],'Box','off','TickDir','out','TickLength',.01*[1 1],'LineWidth',lw,'FontSize',fs)    
set([axh_b; axh_bd; axh_btheta; axh_bphi],'XTickLabel',[])    
set([ph_b; ph_bd; ph_btheta; ph_bphi; ph_omega],'LineWidth',.5,'Color','k')    
set(axh_b,'YLim',max(xps.b/1e9)*[-.05 1.05])    
set(axh_bd,'YLim',[-.55 1.05])    
set(axh_btheta,'YLim',180*[-.05 1.05])    
set(axh_bphi,'YLim',180*[-1.05 1.05])    
set(axh_omega,'YLim',max(xps.rmsomega/2/pi)*[-.05 1.05])    

set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize_p],'PaperSize', papersize_p); 
eval(['print ' voxel_path '/protocol -loose ' pformat])

%%
left_d = .15;
bottom_d = .1;
width_d = .85;
height_d = .8;



nk_v = nk*ones(size(ni_v));
hcontour_v = [];
for npixel = 1:numel(ni_v)
    ni = ni_v(npixel);
    nj = nj_v(npixel);
    nk = nk_v(npixel);
    
    signal = squeeze(I(ni,nj,nk,:));
    m = squeeze(m_bsmerge(ni,nj,nk,:));
    s0 = squeeze(pmaps.s0(ni,nj,nk));
    s_fit = feval([method '_1d_fit2data'], m, xps);
    
    figure(3), clf
    axh = axes('position',[left_d bottom_d width_d height_d]); hold on
    h1 = plot(axh,1:xps.n,signal(ind_sort)/s0,'o');
    h1_fit = plot(axh,1:xps.n,s_fit(ind_sort)/s0,'k.');
    set(h1,'MarkerSize',ms,'LineWidth',.5)
    set(h1,'Color',col_a(npixel,:))
    axis(axh,'tight')
    set(axh,'XLim',xps.n*[-.05 1.05], 'YLim',[-.1 1.1],...
        'Box','off','TickDir','out','TickLength',.01*[1 1],...
    'FontSize',fs,'LineWidth',lw)
    set(axh,'YLim',max([s_fit(:)/s0; signal(:)/s0])*[-.1 1.1])
%     xlabel(axh,'Acq number','FontSize',fs)
%     ylabel(axh,'Signal','FontSize',fs)

    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize_s],'PaperSize', papersize_s); 
    eval(['print ' voxel_path '/signal_' num2str(ni) '_' num2str(nj) '_' num2str(nk) ' -loose ' pformat])
    
%     a = squeeze(I(ni,nj,nk,:));  a = abs(a(:))+eps; a(~isfinite(a)) = 0;
%     clear c
%     c.r = abs(xps.u(:,1)); c.r = c.r(:);
%     c.g = abs(xps.u(:,2)); c.g = c.g(:);
%     c.b = abs(xps.u(:,3)); c.b = c.b(:);
%     c.bright = (abs(bpar-bperp)./max([bpar bperp],[],2)).^2; c.bright = repmat(c.bright(:),[1 3]);
% 
%     figure(3), clf
%     axh = axes('position',[left bottom width height]);
%     ph = scatter3(x,y,z,.2*a,c.bright.*[c.r c.g c.b]);
%     view(30,30)
%     xmin = 0; xmax = 2; ymin = -.5; ymax = 1; zmin = 0; zmax = 65e-3; 
%     axis([xmin xmax ymin ymax zmin zmax])
%     set(ph,'LineWidth',lw_plot)
%     set(axh,'LineWidth',lw_axes,'FontSize',fs_axes)
%     axis(axh,'square')
%     set(axh,'XTick',0:.5:2,'YTick',-.5:.5:1,'ZTick',0:.05:.2)
% 
%     set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
%     eval(['print ' voxel_path '/signal_' num2str(ni) '_' num2str(nj) '_' num2str(nk) ' -loose -dpdf'])


    fh_spec = figure(1); clf
    hax = axes('position',[.25 .1 .65 .65]); hold on
    hax_projx = axes('position',[.25 .75 .65 .25]); hold on
    hax_projy = axes('position',[0 .1 .25 .65]); hold on
%     hold([hax; hax_projx; hax_projy],'on')
    
    
    dtod = dtod_par2dist(squeeze(dpar_4d(ni,nj,nk,:)),squeeze(dperp_4d(ni,nj,nk,:)),squeeze(theta_4d(ni,nj,nk,:)),...
        squeeze(phi_4d(ni,nj,nk,:)),squeeze(d0_4d(ni,nj,nk,:)),squeeze(rpar_4d(ni,nj,nk,:)),squeeze(rperp_4d(ni,nj,nk,:)),...
        squeeze(w_4d(ni,nj,nk,:)));
    for nomega = 1:numel(omega_v)
        [n,dparo,dperpo,theta,phi,w] = dtod_dist2parso(dtod,omega_v(nomega));
        [disoo,danisoo,dratioo,ddeltao,sdanisoo,sddeltao] = dtd_pars2dpars(dparo,dperpo);
%         dist_d.x = log10(disoo);
%         dist_d.y = log10(dratioo);
        dist_d.x = disoo/1e-9;
        dist_d.y = sddeltao;
        dist_d.a = 1000*w/max(w(:));
        dist_d.r = abs(cos(phi).*sin(theta));
        dist_d.g = abs(sin(phi).*sin(theta));
        dist_d.b = abs(cos(theta));
        dist_d.bright = (abs(squeeze(dparo)-squeeze(dperpo))./max([squeeze(dparo) squeeze(dperpo)],[],2)).^2;

        contourpars.Nx = 64; contourpars.Ny = contourpars.Nx; contourpars.Nlevels = 3;
%         axpars.xmin = -11; axpars.xmax = -8; axpars.ymin = -3; axpars.ymax = 3; 
        axpars.xmin = -.2; axpars.xmax = 2.2; axpars.ymin = -.1; axpars.ymax = 1.1; 
        axpars.no_scatter = 1;

        [hax,hscatter,hcontour,hprojx,hprojy] = dist_2d_scattercontourprojplot(hax,hax_projx,hax_projy,dist_d,contourpars,axpars);
%         col = [omega_v(nomega)/max(omega_v) 0 1-omega_v(nomega)/max(omega_v)];
        col = .8*omega_v(nomega)/max(omega_v) + [0 0 0];
        lw_line = 1*lw*(1 - .8*omega_v(nomega)/max(omega_v));

        for nc = 1:numel(hcontour)
            set(hcontour(nc),'Color',col,'LineWidth',lw_line)
        end
        set([hprojx; hprojy],'Color',col,'LineWidth',2*lw_line)
        hcontour_v = [hcontour_v; hcontour];

    end

%     set(hax,'YTick',-4:2:4)
    set(hax,'YTick',0:.5:1,'TickDir','out','TickLength',.02*[1 1],'Box','on')
    set([hax; hax_projx; hax_projy],'FontSize',fs,'LineWidth',lw)
    xlim = get(hax_projx,'XLim'); axis(hax_projx,'tight'); ylim = get(hax_projx,'YLim'); set(hax_projx,'XLim',xlim,'YLim',max(ylim)*[-.1 1.1])
    ylim = get(hax_projy,'YLim'); axis(hax_projy,'tight'); xlim = get(hax_projy,'XLim'); set(hax_projy,'YLim',ylim,'XLim',max(xlim)*[-.1 1.1])

    
    set(fh_spec, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize_d],'PaperSize', papersize_d); 
    eval(['print ' voxel_path '/' method '_' num2str(ni) '_' num2str(nj) '_' num2str(nk) ' -loose ' pformat])
    
    
%     [dxx,dyy,dzz,dxy,dxz,dyz] = dtd_pars2elements(dpar,dperp,theta,phi);
%     [diso,daniso,dratio,ddelta,sdaniso,sddelta] = dtd_pars2dpars(dpar,dperp);
%     
%     figure(4), clf
%     hax = axes('position',[left_d bottom_d width_d height_d]);
% 
%     dist_d.x = log10(diso(ni,nj,nk,:));
%     dist_d.y = log10(dratio(ni,nj,nk,:));
%     dist_d.z = log10(r2(ni,nj,nk,:));
%     dist_d.a = 1*w(ni,nj,nk,:);
%     dist_d.r = abs(cos(phi(ni,nj,nk,:)).*sin(theta(ni,nj,nk,:)));
%     dist_d.g = abs(sin(phi(ni,nj,nk,:)).*sin(theta(ni,nj,nk,:)));
%     dist_d.b = abs(cos(theta(ni,nj,nk,:)));
%     dist_d.bright = (abs(squeeze(dpar(ni,nj,nk,:))-squeeze(dperp(ni,nj,nk,:)))./max([squeeze(dpar(ni,nj,nk,:)) squeeze(dperp(ni,nj,nk,:))],[],2)).^2;
% 
%     contourpars.Nx = 50; contourpars.Ny = contourpars.Nx; contourpars.Nz = contourpars.Nx; contourpars.Nlevels = 3;
%     axpars.xmin = -10; axpars.xmax = -8; axpars.ymin = -3; axpars.ymax = 3; axpars.zmin = 0; axpars.zmax = 2.5; 
% 
%     [hax,~,~] = dist_3d_scattercontourplot(hax,dist_d,contourpars,axpars);
% 
%     set(hax,'LineWidth',lw_axes,'FontSize',fs_axes)
%     set(hax,'XTick',-11:.5:-8,'YTick',-3:1:3,'ZTick',0:.5:2.5)
%     set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
%     eval(['print ' voxel_path '/dtr2d_' num2str(ni) '_' num2str(nj) '_' num2str(nk) ' -loose ' pformat])



%     figure(5), clf
%     hax = axes('position',[left bottom width height]);
% 
%     dist_d.x = log10(dpar(ni,nj,nk,:));
%     dist_d.y = log10(dperp(ni,nj,nk,:));
%     dist_d.z = log10(r2(ni,nj,nk,:));
% 
%     axpars.xmin = -10; axpars.xmax = log10(5e-9); axpars.ymin = axpars.xmin; axpars.ymax = axpars.xmax; axpars.zmin = 0; axpars.zmax = log10(30); 
% 
%     [hax,~,~] = dist_3d_scattercontourplot(hax,dist_d,contourpars,axpars);
% 
%     set(hax,'LineWidth',lw_axes,'FontSize',fs_axes)
%     set(hax,'XTick',-11:.5:-8,'YTick',-11:.5:-8,'ZTick',0:.25:2)
% 
% 
%     set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
%     eval(['print ' paths.voxel_path '/dtr2d_parperp_' num2str(ni) '_' num2str(nj) '_' num2str(nk) ' -loose -dpdf'])
% 
%     figure(6), clf
%     hax = axes('position',[left bottom width height]);
% 
%     dist_d.x = diso(ni,nj,nk,:)/1e-9;
%     dist_d.y = sqdanison(ni,nj,nk,:);
%     dist_d.z = r2(ni,nj,nk,:);
% 
%     axpars.xmin = 0; axpars.xmax = 4; axpars.ymin = 0; axpars.ymax = 1; axpars.zmin = 0; axpars.zmax = 30; 
% 
%     [hax,~,~] = dist_3d_scattercontourplot(hax,dist_d,contourpars,axpars);
% 
%     set(axh,'LineWidth',lw_axes,'FontSize',fs_axes)
%     set(hax,'XTick',0:1:4,'YTick',0:.25:1,'ZTick',0:5:30)
% 
%     set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
%     eval(['print ' paths.voxel_path '/dtr2d_isosqdelta_' num2str(ni) '_' num2str(nj) '_' num2str(nk) ' -loose -dpdf'])

end
%%
delete([hax_projx; hax_projy; hcontour_v])
opt.dtd.bin_disomin = [0 0 1]*1e-9; opt.dtd.bin_disomax = [1 1 5]*1e-9;
opt.dtd.bin_sddeltamin = [.25 0 0]; opt.dtd.bin_sddeltamax = [1 .25 1];    

hold on
Nbin = numel(opt.dtd.bin_disomax);
for nbin = 1:Nbin
    disomin = opt.dtd.bin_disomin(nbin);
    disomax = opt.dtd.bin_disomax(nbin);
    sddeltamin = opt.dtd.bin_sddeltamin(nbin);
    sddeltamax = opt.dtd.bin_sddeltamax(nbin);

%     ph = plot(hax,([disomin disomax disomax disomin disomin]/1e-9),...
%         ([sddeltamin sddeltamin sddeltamax sddeltamax sddeltamin]),...
%         'k-','LineWidth',.5*lw);
    ph = patch(hax,'XData',([disomin disomax disomax disomin disomin]/1e-9),...
        'YData',([sddeltamin sddeltamin sddeltamax sddeltamax sddeltamin]));
    set(ph,'FaceColor',col_a(nbin,:))
end

set(fh_spec, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize_d],'PaperSize', papersize_d); 
eval(['print ' voxel_path '/' method '_bins -loose ' pformat])




return
%%

%Plot bin-resolved global distributions
dist_s.x = linspace(min(r2min),max(r2max),200)';
dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));

figure(7), clf
dist_d.x = mr2_comp3(:);
dist_d.w = mask(:).*f_comp3(:);
dist_d.n = numel(dist_d.x);
dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
plot(dist_s.x,dist_s.w,'-r')
hold on
dist_d.x = mr2_comp2(:);
dist_d.w = mask(:).*f_comp2(:);
dist_d.n = numel(dist_d.x);
dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
plot(dist_s.x,dist_s.w,'-g')
dist_d.x = mr2_comp1(:);
dist_d.w = mask(:).*f_comp1(:);
dist_d.n = numel(dist_d.x);
dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
plot(dist_s.x,dist_s.w,'-b')

dist_s.x = linspace(min(disomin),max(disomax),200)';
dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));

figure(8), clf
dist_d.x = mdiso_comp3(:);
dist_d.w = mask(:).*f_comp3(:);
dist_d.n = numel(dist_d.x);
dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
plot(dist_s.x,dist_s.w,'-r')
hold on
dist_d.x = mdiso_comp2(:);
dist_d.w = mask(:).*f_comp2(:);
dist_d.n = numel(dist_d.x);
dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
plot(dist_s.x,dist_s.w,'-g')
dist_d.x = mdiso_comp1(:);
dist_d.w = mask(:).*f_comp1(:);
dist_d.n = numel(dist_d.x);
dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
plot(dist_s.x,dist_s.w,'-b')

dist_s.x = linspace(min(sqdanisonmin),max(sqdanisonmax),200)';
dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));

figure(9), clf
dist_d.x = msqdanison_comp3(:);
dist_d.w = mask(:).*f_comp3(:);
dist_d.n = numel(dist_d.x);
dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
plot(dist_s.x,dist_s.w,'-r')
hold on
dist_d.x = msqdanison_comp2(:);
dist_d.w = mask(:).*f_comp2(:);
dist_d.n = numel(dist_d.x);
dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
plot(dist_s.x,dist_s.w,'-g')

dist_d.x = msqdanison_comp1(:);
dist_d.w = mask(:).*f_comp1(:);
dist_d.n = numel(dist_d.x);
dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
plot(dist_s.x,dist_s.w,'-b')




