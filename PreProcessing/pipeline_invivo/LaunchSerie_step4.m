%% Launch serie of step 4
clearvars; close all; clc;

addpath('supplementary_functions')
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

%% Define path list
Root_path = 'C:\Users\User\Documents\Data_kuopio\MMD_invivo_final';
path_list = {'ON-101\ses-invivo\20230601_150842_ON_101_3d_1_1\23'};
%   
% 'ON-83/ses-invivo/20230511_090718_ON_83_3d_1_1/13',...
%     'ON-84/ses-invivo/20230511_115103_ON_84_3d_1_1/14',...
%     'ON-85/ses-invivo/20230514_095602_ON_85_3d_1_1/17',...
%     'ON-86/ses-invivo/20230514_123254_ON_86_3d_1_1/28',...
%     'ON-87/ses-invivo/20230515_102145_ON_87_3d_1_1/11',...
%     'ON-88/ses-invivo/20230515_133135_ON_88_3d_run02_1_1/15',...
%     'ON-89/ses-invivo/20230519_093825_ON_89_3d_1_1/20',...
%     'ON-90/ses-invivo/20230519_125904_ON_90_3d_1_1/12',...
%     'ON-91/ses-invivo/20230520_091240_ON_91_3d_1_1/13',...
%     'ON-92/ses-invivo/20230520_114937_ON_92_3d_1_1/13',...
%     'ON-93/ses-invivo/20230520_143622_ON_93_3d_1_1/16',...
%     'ON-94/ses-invivo/20230525_101903_ON_94_3d_1_1/14',...
%     'ON-95/ses-invivo/20230525_132517_ON_95_3d_1_1/26',...
%     'ON-96/ses-invivo/20230527_083012_ON_96_3d_1_1/12',...
%     'ON-97/ses-invivo/20230527_110209_ON_97_3d_1_1/19',...
%     'ON-98/ses-invivo/20230527_140830_ON_98_3d_1_1/16'};
         
%% Pre-processing loop
for ind = 1:size(path_list,2)
    disp(['Started exp ' num2str(ind)])
    %% set current data path
    data_path = path_list{ind};
    data_path = [Root_path filesep data_path filesep 'pdata_mdd'];
    data_path = split(data_path,'/');
    data_path = join(data_path(1:end,1),filesep,1); data_path = data_path{1};
    function_step4_plots_batch_invivo_rat(data_path)
    disp(['Finsished exp ' num2str(ind)])
end
