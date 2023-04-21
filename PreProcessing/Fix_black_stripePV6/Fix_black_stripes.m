%% Fix Reco Overflow error
clearvars; close all; clc;

%% Define paths
pathData = 'C:\Users\User\Mon_Drive\Matlab\ProcessingPV360\data\ON-79-mtbi-pilot-run02\9\pdata\1\';

%% Save original visu_par file
if ~isfile([pathData 'saved_visu_pars'])
   copyfile([pathData 'visu_pars'],[pathData 'saved_visu_pars'],'f');
else
   copyfile([pathData 'saved_visu_pars'],[pathData 'visu_pars'],'f');  
end

%% open VisuCoreDataSlope and get value
[visu_pars,~]=readBrukerParamFile([pathData 'visu_pars']);

% % Display
figure(1)
plot(visu_pars.VisuCoreDataSlope)
title('I should look like an horizontal line with few negative spikes')

VisuCoreDataSlopeVal = median(visu_pars.VisuCoreDataSlope);
New_VisuCoreDataSlope = repmat(VisuCoreDataSlopeVal,1,size(visu_pars.VisuCoreDataSlope,2));
%% Fix value
modify_visu_pars(pathData,'VisuCoreDataSlope',New_VisuCoreDataSlope);

disp('VisuCoreDataSlope parameter fixed')


