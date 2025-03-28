function [] = my_function_step2_pv2nii_xps(data_path)
% Convert Bruker ParaVision data to nii and xps

% Image recon parameters
rps.denoising       = 0; %denoising
rps.maxomega        = 300*2*pi;
rps.nomega          = 30;

%-------------------------------------------------
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
procno =1;
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
    if strcmp(Method_name,'<user:my_dor_r12_epi>' )==1 || strcmp(Method_name,'<user:my_dor_r12_epi_cs>')==1
        mdm_bruker_my_DOR_R12_EPI_PV360_acqus2xps(data_path, xps_temp_fn, rps);
    end
    if strcmp(Method_name,'<user:my_dor_r12_epi_shift>')==1 
        mdm_bruker_my_DOR_R12_EPI_Shift_PV360_acqus2xps(data_path, xps_temp_fn, rps);
    end
    if strcmp(Method_name,'<user:my_dor_r12_epi_comp>')==1
        mdm_bruker_my_DOR_R12_EPI_comp_PV360_acqus2xps(data_path, xps_temp_fn, rps);
    end
    if strcmp(Method_name,'<user:my_dor_r12_epi_now>')==1
        mdm_bruker_my_DOR_R12_EPI_NOW_PV360_acqus2xps(data_path, xps_temp_fn, rps);
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

close all; clc;


