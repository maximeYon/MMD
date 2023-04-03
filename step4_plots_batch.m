%Plot parameters maps and global stats
clearvars; close all; clc;

home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

% Define full paths to the nifti files to be analyzed 
% data_path = '20211019_110649_test_sample_Hong_test_sample_Hong_1_1\127';
% data_path = '20220204_tumor3\60';
% data_path = '20220202_RatBrain\22';
% data_path = '20220208_HongPhantom3\9';
% data_path = 'Original_Rat_Brain\14'; 
% data_path = 'Cerebelum\18'; 
% data_path = '20220407_Mouse_heart_PBS\65';
% data_path = '20220411_Mouse_Brain_1\32'; 
% data_path = '20220420_Mouse_brain2\15'; 
% data_path = '20220524_mouse_brain_fin2\53';
% data_path = '20221103_mouse_brain_fin7\22'; 
% data_path = 'omar_hippocampus_sample\10';
% data_path = '20221211_mouse_heart\22';
% data_path = '20221206_kenneth_data\30'; 
% data_path = 'human_brain\1'; 
% data_path = '20230118_hippo3D\17';
data_path = 'human_sample_fresh1\34';

data_path = split(data_path,'\'); data_path = join(data_path,filesep,1); data_path = data_path{1};
nii_fns{1} = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path filesep 'pdata_mdd' filesep 'nii_xps' filesep 'data.nii.gz'];

%-----------------------------------------
% Make paths to nii_xps, boostraps, and maps folders

% method = 'dtd';
% method = 'dtod'; omega_v = [53 160]*2*pi;
% method = 'dtr1r2d';
%method = 'dtor1r2d'; omega_v = [12.7 387]*2*pi; % 10% -> 35 Hz
% method = 'dtor1r2d'; omega_v = [10 70]*2*pi; % in vivo at 7T
% method = 'dtor1r2d'; omega_v = [100 387]*2*pi; % protocol 303 [12.7 387]
method = 'dtor1r2d'; omega_v = [23 332]*2*pi; % protocol 303 4 shapes, mean = 173
% method = 'dtor1r2d'; omega_v = [9.04 11]*2*pi; % humain brain: 10 % : 5.15 90% 11 mean 9.04

pmaps_paths = cell(0,0);
bs_paths = cell(0,0);
for ndata = 1:numel(nii_fns)
    %pdata_path = fullfile(fileparts(fileparts(nii_fns{ndata})),'opt1');
    pdata_path = fileparts(fileparts(nii_fns{ndata}));
    pmaps_paths{1+numel(pmaps_paths)} = fullfile(pdata_path,'pmaps');
    bs_paths{1+numel(bs_paths)} = fullfile(pdata_path,method,'bootstraps');        
end

% Define bins for ex vivo rat brain
disomin = [0 0 1]*1e-9; disomax = [1 1 5]*1e-9;
sddeltamin = [.25 0 0]; sddeltamax = [1 .25 1];
dratiomin = ones(size(disomin))*eps; dratiomax = ones(size(disomin))/eps;
r1min = 0*ones(size(disomin)); r1max = 21*ones(size(disomin));
r2min = 0*ones(size(disomin)); r2max = 201*ones(size(disomin));

% Define bins for in vivo human brain
% disomin = [0 0 2.5]*1e-9; disomax = [1 1 5]*1e-9;
% sddeltamin = [.25 0 0]; sddeltamax = [1 .25 1];
% dratiomin = ones(size(disomin))*eps; dratiomax = ones(size(disomin))/eps;
% r1min = 0*ones(size(disomin)); r1max = 21*ones(size(disomin));
% r2min = 0*ones(size(disomin)); r2max = 201*ones(size(disomin));

% Define color limits for parameter maps
clim.s0 = .9*[0 1]; % Multiplied with s0 below
clim.mdiso = 1e-9*[0 1]; %1e-9*[0 1]
clim.msddelta = 1*[0 1];
clim.mr1 = 1*[0 1]; %1
clim.mr2 = 50*[0 1];
clim.vdiso = .8e-18*[0 1];
clim.vsddelta = .20*[0 1];
clim.vr1 = 0.2*[0 1];
clim.vr2 = 500*[0 1];
clim.cvdisosddelta = 0.8*sqrt(max(clim.vdiso)*max(clim.vsddelta))*[-1 1]; %clim.cvdisosddelta = .2e-9*[-1 1];
clim.cvdisor1 = 0.4*sqrt(max(clim.vdiso)*max(clim.vr1))*[-1 1];
clim.cvdisor2 = 0.5*sqrt(max(clim.vdiso)*max(clim.vr2))*[-1 1];
clim.cvsddeltar1 = 1*sqrt(max(clim.vsddelta)*max(clim.vr1))*[-1 1];
clim.cvsddeltar2 = 1*sqrt(max(clim.vsddelta)*max(clim.vr2))*[-1 1];
clim.cvr1r2 = 0.5*sqrt(max(clim.vr1)*max(clim.vr2))*[-1 1];
clim.mask_threshold = eps; %clim.mask_threshold = 0.01;

clim.dmdisodnu = 2e-3*max(clim.mdiso)*[-1 1];
% clim.dmdisodnu = 1e-3*max(clim.mdiso)*[0 1];
clim.dmsddeltadnu = 1e-3*max(clim.msddelta)*[-1 1];
clim.dvdisodnu = 1e-2*max(clim.vdiso)*[-1 1];
clim.dvsddeltadnu = 1e-3*max(clim.vsddelta)*[-1 1];
clim.dcvdisosddeltadnu = 1e-3*max(clim.cvdisosddelta)*[-1 1];

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









