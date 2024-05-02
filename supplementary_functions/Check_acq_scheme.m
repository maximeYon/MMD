clearvars;
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

% data_path = '20211220_tumor2\276';%'20211220_tumor2\27' 
% data_path = '20220125_112215_Test_half_brain_dtor1r2d_maxime_01_Test_hal_1_1\43'; 
% data_path = '20220620_130815_qMAS_testing\132\';
data_path ='clay_1\13'; 
data_path = split(data_path,'\'); data_path = join(data_path,filesep,1); data_path = data_path{1};
data_path = [home_path filesep 'ProcessingPV360' filesep 'data' filesep data_path];

s.nii_fn = [data_path '\pdata_mdd\nii_xps\data.nii.gz'];

% Add sorting by Bvalue & Bdelta option
my_mdm_niixps2pdf(s,'Yes') % Yes
% mdm_niixps2pdf(s)