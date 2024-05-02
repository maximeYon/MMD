%% Open PRESS spectra
clearvars; close all; clc;
file_path = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\20220620_130815_qMAS_testing\136\pdata_mdd\nii_xps';

% Load data
opt = mdm_mrtrix_opt;
[data,h] = mdm_nii_read([file_path filesep 'data.nii.gz']);
data = squeeze(data);

% Load xps
load([file_path filesep 'data_xps.mat'])

% display
figure(1)
mesh(data)

% Integrate resonances
Nres = 3;
coord = [1365,1455;1692,1775;1945,2142];
for ind = 1:Nres
    S_int{ind} = sum(data(coord(ind,1):coord(ind,2),:),1);
end

% display decay curves
figure(2)
hold on
for ind = 1:Nres
    plot(xps.b/10^9,S_int{ind}/max(S_int{ind}),'*')
end







