%% Verify the maximum gradient amplitude in R P S in DTD dor experiment
clearvars; close all; clc;

%Load list
% path_list = '/opt/PV-360.1.1/prog/curdir/maxime/ParaVision/exp/lists/DTD/DOR_R1R2_1404';
path_list = 'C:\Users\User\Mon_Drive\Matlab\MATLABfunctions\md-dmri\acq\bruker\Paravision360\DOR_list\DOR_R1R2_303_4';
LowLevelDTDr1r2list = importfileDTDlist(path_list, [3 inf]);
LowLevelDTDr1r2list = table2array(LowLevelDTDr1r2list(1:end-1,:));

%Load Frequency Shapes
% path_shapes = '/opt/PV-360.1.1/prog/curdir/maxime/ParaVision/exp/lists/gp/';
path_shapes = 'C:\Users\User\Mon_Drive\Matlab\MATLABfunctions\md-dmri\acq\bruker\Paravision360\DOR_shapes\';

Nshapes = 6;
DORShapStruct = struct;
for ind = 1:Nshapes
DORShapStruct.(['shape' num2str(ind)]) = table2array(importfile_DOR_Shape([path_shapes 'DOR_freq' num2str(ind) '_1000pts'],[22, Inf]));
DORShapStruct.(['shape' num2str(ind)]) = DORShapStruct.(['shape' num2str(ind)])(1:end-1,:);
end
clearvars path_list path_shapes ind;

%% Compute gradients
Gx = zeros(size(DORShapStruct.shape1,1),size(LowLevelDTDr1r2list,1));
Gy = zeros(size(DORShapStruct.shape1,1),size(LowLevelDTDr1r2list,1));
Gz = zeros(size(DORShapStruct.shape1,1),size(LowLevelDTDr1r2list,1));
for ind = 1:size(LowLevelDTDr1r2list,1)
    shape_number = LowLevelDTDr1r2list(ind,2)+1;
    if shape_number>0
    Gx(:,ind) = LowLevelDTDr1r2list(ind,4)*DORShapStruct.(['shape' num2str(shape_number)])(:,1) +...
            LowLevelDTDr1r2list(ind,5)*DORShapStruct.(['shape' num2str(shape_number)])(:,2) +...
            LowLevelDTDr1r2list(ind,6)*DORShapStruct.(['shape' num2str(shape_number)])(:,3);
        
     Gy(:,ind) = LowLevelDTDr1r2list(ind,7)*DORShapStruct.(['shape' num2str(shape_number)])(:,1) +...
              LowLevelDTDr1r2list(ind,8)*DORShapStruct.(['shape' num2str(shape_number)])(:,2) +...
              LowLevelDTDr1r2list(ind,9)*DORShapStruct.(['shape' num2str(shape_number)])(:,3);
          
     Gz(:,ind) = LowLevelDTDr1r2list(ind,10)*DORShapStruct.(['shape' num2str(shape_number)])(:,1) +...
              LowLevelDTDr1r2list(ind,11)*DORShapStruct.(['shape' num2str(shape_number)])(:,2) +...
              LowLevelDTDr1r2list(ind,12)*DORShapStruct.(['shape' num2str(shape_number)])(:,3);  
    end
end

%% Display
figure()
subplot(3,1,1)
plot(Gx)
subplot(3,1,2)
plot(Gy)
subplot(3,1,3)
plot(Gz)

max([Gx(:) Gy(:) Gz(:)])

figure()
Gtot = Gx.^2+Gy.^2+Gz.^2;
plot(Gtot)
max(Gtot(:))
% aa = (1/sqrt(3)).^2+(1/sqrt(3)).^2+(1/sqrt(3)).^2;




