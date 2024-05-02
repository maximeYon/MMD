% Convert Bruker ParaVision data to nii and xps
clearvars;close all; clc;

% Add function to path
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

% If spectroscopy
addpath(genpath([home_path filesep 'ProcessingPV360' filesep 'spectroscopy_diff' filesep 'spectroscopy_functions' ]));

% Image recon parameters
rps.denoising       = 0; %denoising
rps.maxomega        = 600*2*pi; % 300 for dtd ; 1200 for the reste
rps.nomega          = 50; %50 for DOR

%-------------------------------------------------
% Define path to the experiment to be analyzed 
% data_path = '20220204_tumor3\60'; 
% data_path = '20220208_HongPhantom3\14';
% data_path = '20220202_RatBrain\22'; 
% data_path = 'fresh_Brain\8'; 
%data_path = 'Cerebelum\18'; 
% data_path = '20220406_Hong_phamtom5\3'; 
% data_path = '20220407_Mouse_heart_PBS\65'; 
%data_path = '20220411_Mouse_Brain_1\24'; 
%data_path = '20220411_Mouse_Brain_1\32'; 
% data_path = '20220420_Mouse_brain2\15'; 
% data_path = '20220421_Tumor5\9'; 
% data_path = '20220411_Mouse_Brain_1\41'; 
% data_path = '20220524_mouse_brain_fin2\73'; 
% data_path = '20220722_mouse_brain_fin4\9'; 
% data_path = 'spectro\37'; 
% data_path = '20221103_mouse_brain_fin7\22'; 
% data_path = 'omar_hippocampus_sample\9'; 
% data_path = '20221206_kenneth_data\30'; 
% data_path = '\20221209_Hong_phantom_spectro\7';
% data_path = '20221211_mouse_heart\14';
% data_path = 'clay_1\13';
% data_path = '20230118_hippo3D\24';
% data_path = '20230206_3D_hippo\40';
% data_path = '3D_human1\46';
data_path = 'invivo_rat4\27';


procno =1;

data_path = split(data_path,'\'); data_path = join(data_path,filesep,1); data_path = data_path{1};
data_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path];

% Create Out Path
out_path = fullfile(data_path, 'pdata_mdd','nii_xps');
msf_mkdir(out_path);

%% Perform conversion
% get opt parameters
opt = mdm_opt();

% Store computed parameter in tmp files
nii_temp_fn = fullfile(out_path, ['data_temp' opt.nii_ext]);
xps_temp_fn = fullfile(out_path, 'data_temp_xps.mat');
%% retrieve method name
Method_name = ReadPV360Param([data_path filesep], 'Method');

%% retrieve Paravision version
PV_version = ReadPV360Param([data_path filesep], 'ACQ_sw_version');

%% Convert data, perform denoising if requiered
mdm_bruker_PV360_2dseq2nii(data_path, nii_temp_fn, rps,procno);

%% Compute and store the b-values
% if contains(Method_name,'diff_wave_sym_msme')==1
%      mdm_bruker_diff_wave_sym_msme_acqus2xps(data_path, xps_temp_fn, rps);
% end
if contains(PV_version,'pv 6.')==1
    if strcmp(Method_name,'<user:my_dor_r12_rare>')==1 || strcmp(Method_name,'<user:my_dor_r12_rare_vfa>')==1 || strcmp(Method_name,'<user:my_dor_r12_rare_cs>')==1 || strcmp(Method_name,'<user:my_dor_r12_rare_full>')==1 || strcmp(Method_name,'<user:my_dorr12_rare_final>')==1
        mdm_bruker_my_DOR_R12_RARE_PV6_acqus2xps(data_path, xps_temp_fn, rps);
    end
    if strcmp(Method_name,'<user:my_dor_r12_msme>')==1
        mdm_bruker_my_DOR_R12_MSME_PV6_acqus2xps(data_path, xps_temp_fn, rps);
    end  
    if strcmp(Method_name,'<user:my_dor12_press>')==1
        mdm_bruker_my_DOR_R12_PRESS_PV6_acqus2xps(data_path, xps_temp_fn, rps);
    end  
    if strcmp(Method_name,'<user:my_dor12_press_e2>')==1
        mdm_bruker_my_DOR_R12_PRESS_e2_PV6_acqus2xps(data_path, xps_temp_fn, rps);
    end  
    if strcmp(Method_name,'<user:my_dor12_epi>')==1 || strcmp(Method_name,'<user:my_dor12_epi_tup>')==1 || strcmp(Method_name,'<user:my_dor12_epi_multsl>')==1 
        mdm_bruker_my_DOR_R12_EPI_PV6_acqus2xps(data_path, xps_temp_fn, rps);
    end  
    if strcmp(Method_name,'<user:my_dor12_spiral>')==1
        mdm_bruker_my_DOR_R12_SPIRAL_PV6_acqus2xps(data_path, xps_temp_fn, rps);
    end  
end

if contains(PV_version,'pv-360')==1
    if strcmp(Method_name,'<user:my_dor_r12_segflash>')==1 
        mdm_bruker_my_DOR_R12_SegFLASH_PV360_acqus2xps(data_path, xps_temp_fn, rps);
    end
    if strcmp(Method_name,'<user:my_dor_r12_rare>')==1 || strcmp(Method_name,'<user:my_dor_r12_rare_vfa>')==1 || strcmp(Method_name,'<user:my_dor_r12_rare_cs>')==1 || strcmp(Method_name,'<user:my_dor_r12_rare_full>')==1 || strcmp(Method_name,'<user:my_dor12_rare_final>')==1 
        mdm_bruker_my_DOR_R12_RARE_PV360_acqus2xps(data_path, xps_temp_fn, rps);
    end
    if strcmp(Method_name,'<user:my_dor_r12_epi>')==1
        mdm_bruker_my_DOR_R12_EPI_PV360_acqus2xps(data_path, xps_temp_fn, rps);
    end
    if strcmp(Method_name,'<user:my_dor_r12_msme>')==1
        mdm_bruker_my_DOR_R12_MSME_PV360_acqus2xps(data_path, xps_temp_fn, rps);
    end
    if strcmp(Method_name,'<user:my_dor_r12_msme_opt>')==1
        mdm_bruker_my_DOR_R12_MSME_opt_PV360_acqus2xps(data_path, xps_temp_fn, rps);
    end
    if strcmp(Method_name,'<user:my_dor12_press>')==1
        mdm_bruker_my_DOR12_PRESS_PV360_acqus2xps(data_path, xps_temp_fn, rps);
    end  
end

%% Move tmpfile to permanent files
% RECO_image_type=ReadPV360Param(data_path_pv, 'RECO_image_type', [data_path_pv 'pdata' filesep num2str(procno) filesep]);
% if strfind(RECO_image_type,'complex')
% end

s.nii_fn = nii_temp_fn;
load(xps_temp_fn);
s.xps = xps;

nii_fn = fullfile(out_path, ['data' opt.nii_ext]);
xps_fn = fullfile(out_path, 'data_xps.mat');

movefile(xps_temp_fn,xps_fn,'f');

% RECO_image_type=ReadPV360Param([data_path filesep], 'RECO_image_type', [data_path filesep 'pdata' filesep num2str(procno) filesep]);
% if strfind(RECO_image_type,'complex')
if exist([nii_temp_fn(1:end-7) '_real' nii_temp_fn(end-6:end)],'file') && exist([nii_temp_fn(1:end-7) '_imag' nii_temp_fn(end-6:end)],'file')
    movefile([nii_temp_fn(1:end-7) '_real' nii_temp_fn(end-6:end)],[nii_fn(1:end-7) '_real' nii_fn(end-6:end)],'f');
    movefile([nii_temp_fn(1:end-7) '_imag' nii_temp_fn(end-6:end)],[nii_fn(1:end-7) '_imag' nii_fn(end-6:end)],'f');
else
    movefile(nii_temp_fn,nii_fn,'f');
end




