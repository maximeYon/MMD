%% Verify the maximum gradient amplitude in R P S in DTD dor experiment
clearvars; close all; clc;

PVM_RiseTime = 0.152; % ms

%Load list
% path_list = '/opt/PV-360.1.1/prog/curdir/maxime/ParaVision/exp/lists/DTD/DOR_R1R2_1404';
path_list = 'C:\Users\User\Mon_Drive\Matlab\MATLABfunctions\md-dmri\acq\bruker\Paravision360\4_7T_Rennes_BGA12S\DOR_list\DOR_R1R2_389_2';
LowLevelDTDr1r2list = importfileDTDlist(path_list, [3 inf]);
LowLevelDTDr1r2list = table2array(LowLevelDTDr1r2list(1:end-1,:));
shapes_durations = LowLevelDTDr1r2list(:,3);

%Load Frequency Shapes
% path_shapes = '/opt/PV-360.1.1/prog/curdir/maxime/ParaVision/exp/lists/gp/';
path_shapes = 'C:\Users\User\Mon_Drive\Matlab\MATLABfunctions\md-dmri\acq\bruker\Paravision360\4_7T_Rennes_BGA12S\DOR_shapes\';

Nshapes = 4;
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

%% Display Shapes intensities
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

%% Display Shapes derivative and compare with slew rate
slew_rate = 1/PVM_RiseTime; % PVM_RiseTime
dGx = [zeros(1,size(Gx,2));diff(Gx)];
dGy = [zeros(1,size(Gy,2));diff(Gy)];
dGz = [zeros(1,size(Gy,2));diff(Gy)];

figure()
for ind_exp = 1:size(Gx,2)
    shape_time = 0:(shapes_durations(ind_exp)/999):shapes_durations(ind_exp);
    dt = shape_time(1,2)-shape_time(1,1);
    subplot(3,1,1)
    hold on
    plot(shape_time,dGx(:,ind_exp)./dt)
    subplot(3,1,2)
    hold on
    plot(shape_time,dGy(:,ind_exp)./dt)
    subplot(3,1,3)
    hold on
    plot(shape_time,dGz(:,ind_exp)./dt)
    
    if ind_exp ==size(Gx,2) 
    subplot(3,1,1)
    hold on
    plot([0 max(shapes_durations)],[slew_rate slew_rate],'r', 'LineWidth',2)
    subplot(3,1,2)
    hold on
    plot([0 max(shapes_durations)],[slew_rate slew_rate],'r', 'LineWidth',2)
    subplot(3,1,3)
    hold on
    plot([0 max(shapes_durations)],[slew_rate slew_rate],'r', 'LineWidth',2)
    end
end

max(max(abs(dGx./(shapes_durations'/1000))))
max(max(abs(dGy./(shapes_durations'/1000))))
max(max(abs(dGz./(shapes_durations'/1000))))

slew_rate



