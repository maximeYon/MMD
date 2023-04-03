%% De Gibbs NIFTI dataset via MRtrix
clearvars;

home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

% Define full paths to the nifti files to be analyzed 
data_path = 'invivoRat3/21'; 
mdd_path = 'pdata_mdd';

data_path = split(data_path,'\'); data_path = join(data_path,filesep,1); data_path = data_path{1};
nii_fns{1} = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path filesep mdd_path filesep 'nii_xps' filesep 'denoised' filesep 'data.nii'];
mask_fns{1} = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path filesep mdd_path filesep 'nii_xps' filesep 'data_mask.nii.gz'];
nii_fns_out{1} = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path filesep mdd_path filesep 'nii_xps' filesep 'data.nii.gz'];

%% Save original data
% nii_folder = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path filesep 'pdata_mdd' filesep 'nii_xps'];
% if ~isfolder([nii_folder filesep 'orig'])
%    mkdir([nii_folder filesep 'orig']);
%    copyfile([nii_folder filesep 'data.nii.gz'],[nii_folder filesep 'orig' filesep 'data.nii.gz'],'f');
%    copyfile([nii_folder filesep 'data_xps.mat'],[nii_folder filesep 'orig' filesep 'data_xps.mat'],'f');
% end

%% perform MRtrix de-Gibbs
opt = mdm_mrtrix_opt;
[status, result] = my_mdm_mrtrix_degibbs(nii_fns{1}, nii_fns_out{1}, opt);
