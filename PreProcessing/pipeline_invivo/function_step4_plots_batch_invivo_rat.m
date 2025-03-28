function [] = function_step4_plots_batch_invivo_rat(data_path)%Plot parameters maps and global stats

nii_fns{1} = [data_path filesep 'nii_xps' filesep 'data.nii.gz'];

%-----------------------------------------
% Make paths to nii_xps, boostraps, and maps folders
 method = 'dtor1r2d'; omega_v = [15 80]*2*pi; % in vivo at 7T


pmaps_paths = cell(0,0);
bs_paths = cell(0,0);
for ndata = 1:numel(nii_fns)
    %pdata_path = fullfile(fileparts(fileparts(nii_fns{ndata})),'opt1');
    pdata_path = fileparts(fileparts(nii_fns{ndata}));
    pmaps_paths{1+numel(pmaps_paths)} = fullfile(pdata_path,'pmaps');
    bs_paths{1+numel(bs_paths)} = fullfile(pdata_path,method,'bootstraps');        
end


% Define bins for in vivo human brain
disomin = [0 0 2.5]*1e-9; disomax = [1 1 5]*1e-9;
sddeltamin = [.25 0 0]; sddeltamax = [1 .25 1];
dratiomin = ones(size(disomin))*eps; dratiomax = ones(size(disomin))/eps;
r1min = 0*ones(size(disomin)); r1max = 21*ones(size(disomin));
r2min = 0*ones(size(disomin)); r2max = 201*ones(size(disomin));

% Define color limits for parameter maps
clim.s0 = .9*[0 1]; % Multiplied with s0 below
clim.mdiso = 2.5e-9*[0 1]; %1e-9*[0 1]
clim.msddelta = 1*[0 1];
clim.mr1 = 1*[0 1]; %1
clim.mr2 = 40*[0 1];
clim.vdiso = .8e-18*[0 1];
clim.vsddelta = .20*[0 1];
clim.vr1 = 0.2*[0 1];
clim.vr2 = 300*[0 1];
clim.cvdisosddelta = 0.8*sqrt(max(clim.vdiso)*max(clim.vsddelta))*[-1 1]; %clim.cvdisosddelta = .2e-9*[-1 1];
clim.cvdisor1 = 0.4*sqrt(max(clim.vdiso)*max(clim.vr1))*[-1 1];
clim.cvdisor2 = 0.5*sqrt(max(clim.vdiso)*max(clim.vr2))*[-1 1];
clim.cvsddeltar1 = 1*sqrt(max(clim.vsddelta)*max(clim.vr1))*[-1 1];
clim.cvsddeltar2 = 1*sqrt(max(clim.vsddelta)*max(clim.vr2))*[-1 1];
clim.cvr1r2 = 0.5*sqrt(max(clim.vr1)*max(clim.vr2))*[-1 1];
clim.mask_threshold = eps; %clim.mask_threshold = 0.01;


clim.dmdisodnu = 1e-3*max(clim.mdiso)*[-1 1];
clim.dmsddeltadnu = 1e-3*max(clim.msddelta)*[-1 1];
clim.dvdisodnu = 0.3e-2*max(clim.vdiso)*[-1 1];
clim.dvsddeltadnu = 0.3e-3*max(clim.vsddelta)*[-1 1];
clim.dcvdisosddeltadnu = 0.3e-3*max(clim.cvdisosddelta)*[-1 1];

%------------------------------

% Prepare options
opt = mdm_opt();
opt = feval([method '_opt'],opt);
opt.(method).bin_disomin = disomin; opt.(method).bin_disomax = disomax;
opt.(method).bin_dratiomin = dratiomin; opt.(method).bin_dratiomax = dratiomax;
opt.(method).bin_sddeltamin = sddeltamin; opt.(method).bin_sddeltamax = sddeltamax;
if strcmp(method,'dtr2d')
    opt.(method).bin_r2min = r2min; opt.(method).bin_r2max = r2max;
elseif strcmp(method,'dtr1d')
    opt.(method).bin_r1min = r1min; opt.(method).bin_r1max = r1max;
elseif strcmp(method,'dtr1r2d')
    opt.(method).bin_r1min = r1min; opt.(method).bin_r1max = r1max;
    opt.(method).bin_r2min = r2min; opt.(method).bin_r2max = r2max;
elseif strcmp(method,'dtod')
    opt.(method).maps_omega = omega_v;
elseif strcmp(method,'dtor1r2d')
    opt.dtod.maps_omega = omega_v;
    opt.(method).maps_omega = omega_v;
    opt.(method).bin_r1min = r1min; opt.(method).bin_r1max = r1max;
    opt.(method).bin_r2min = r2min; opt.(method).bin_r2max = r2max;
end
opt.mplot.zoomx = 80/100;
opt.mplot.zoomy = 80/100;

Ndata = numel(bs_paths);
tic

% Loop over datasets
for ndata = 1:Ndata
    
    bs_dps = mdm_dps_collectbs(method, bs_paths{ndata}, opt);
        
    if ~all(cellfun('isempty',bs_dps))
        median_dps = mdm_dps_median(bs_dps);
        clear bs_dps
%%  
        mplot_technicolor_nii(method, median_dps, pmaps_paths{ndata}, clim, opt)
        mplot_technicolor_slicemontage2(method, median_dps, fullfile(fileparts(pmaps_paths{ndata}),'slicemontage'), clim, opt)
        mplot_globalstats(method, median_dps, fullfile(fileparts(pmaps_paths{ndata}),[method '_globalstats']), clim)
    end
   
end
toc

% %% Gif creation
% fraction = cat(4,median_dps.bin{1, 1}.f,median_dps.bin{1, 2}.f,median_dps.bin{1, 3}.f);
% FoV = [16 9.5];
% x_im = 0:FoV(1,1)/(size(median_dps.msddelta,1)-1):FoV(1,1);
% y_im =flip(0:FoV(1,2)/(size(median_dps.msddelta,2)-1):FoV(1,2));
% filename = 'human_fresh1_200.gif';
% 
% figure(100)
% set(gcf,'color','w', 'Position', [744 95.4  449  954.6]);
% for Nslice = 1:64
% % Nslice = 32;
% a1 = subplot(3,1,1);
% imagesc(x_im,y_im,flip(flip(squeeze(fraction(:,:,Nslice,:)),2),1));
% axis image
% set(gca, 'Ydir', 'Normal')
% xticks([0 2 4 6 8 10 12 14 16]);yticks([0 2 4 6 8])
% set(gca, 'box','off');
% xlabel('FoV (mm)'); ylabel('FoV (mm)');
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold')
% title('Fractions')
% 
% a2 = subplot(3,1,2);
% imagesc(x_im,y_im,flip(flip(median_dps.msddelta(:,:,Nslice),2),1),[0 1]);
% axis image
% colormap(a2,gray)
% set(gca, 'Ydir', 'Normal')
% xticks([0 2 4 6 8 10 12 14 16]);yticks([0 2 4 6 8])
% set(gca, 'box','off');
% xlabel('FoV (mm)'); ylabel('FoV (mm)');
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold')
% title('E[D_{\Delta} ^{2}]')
% 
% a3 = subplot(3,1,3);
% imagesc(x_im,y_im,flip(flip(median_dps.dmdisodnu(:,:,Nslice),2),1),[0 2e-12]);
% axis image
% colormap(a3,hot(128))
% set(gca, 'Ydir', 'Normal')
% xticks([0 2 4 6 8 10 12 14 16]);yticks([0 2 4 6 8])
% set(gca, 'box','off');
% xlabel('FoV (mm)'); ylabel('FoV (mm)');
% set(gca,'FontSize',12)
% set(gca,'FontWeight','bold')
% title('\Delta_{\omega/2\pi}E[D_{iso}]/10^{-9}')
% 
% drawnow
% pause(0.1)
%       frame = getframe(gcf);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if Nslice == 1
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%       end
% 
% end
% 









