% Acquisition protocol for R1 R2 DOR encoding.
%
% Used in Topgaard, Phys. Chem. Chem. Phys. 18, 8545 (2016).
% http://dx.doi.org/10.1039/c5cp07251d

clearvars;close all;
xlsx_name = 'DOR_protocol_invivo_1-2shapes_389_shortT1'; % 'List_testBvalue'
% addpath(genpath('/home/maxime/Matlab/MATLABfunctions'));

% Define path for output running and repulsion folders
run_path = fileparts(mfilename('fullpath'));
% out_path = '/opt/PV-360.1.1/prog/curdir/maxime/ParaVision/exp/lists/DTD/';
out_path = 'C:\Users\User\Mon_Drive\Matlab\MATLABfunctions\md-dmri\acq\bruker\Paravision360\4_7T_Rennes_BGA12S\DOR_list\';
framework_path = fileparts(fileparts(fileparts(fileparts(run_path)))); % folder dependent
pa_path = fullfile(framework_path,'tools','uvec','repulsion_angles');

%% Define parameters
randomised =1; % randomized option
multislice = 0; % if multisclice =1, sorting is different to have a short TR as first experiment

% Gradients and shapes constant
% PVM_GradCalConst = 128803; % Kuopio 500 Micro MRI
% PVM_GradCalConst = 32597.4; % Kuopio's 7T
PVM_GradCalConst = 28218.7; % Rennes 4.7 T BGA 12S
gamma1H = 26.75e7;
param.Gmax= PVM_GradCalConst*2*pi/gamma1H*1e3;
param.EffiFactor =[0.0123113312115505,0.00360553646759647,0.00135234506482497,0.000702017885394340]; % 1000 points

% tau = 4;
% bval = (2*param.EffiFactor*(1.*param.Gmax)^2*(tau*10^-3)^3*gamma1H^2)*10^-9;
param.Min_Shape_Dur = 4; % minimal duration of shapes in ms

%% read xlsx file
% Import the data
CellParam = readcell([xlsx_name '.xlsx']);
CellParam = CellParam(1:9,2:end);
for c = 1:size(CellParam,2)
    for l = 1:size(CellParam,1)
        CellParam{l,c} = split(replace(string(CellParam{l,c}),",","."),"_");
        CellParam{l,c} = str2double(CellParam{l,c})';
    end
end

%% Get max shape order for the name
max_shape_order_name = 0;
for c = 1:size(CellParam,2)
    max_shape_order_name = max([max_shape_order_name max(CellParam{6,c})]);
end


%% Create array
for c = 1:size(CellParam,2)
    %% column normalization
    max_shape_order = max(CellParam{6,c});
    max_shape_order = max([max_shape_order 1]); % avoid error for no diff
    param.EffiFactorNorm = param.EffiFactor(max_shape_order);
    %% end of column normalization
    
    [Par_Array_supp] = create_array(pa_path,CellParam{1,c}, CellParam{2,c}, CellParam{3,c},CellParam{4,c},CellParam{5,c},CellParam{6,c},CellParam{7,c},CellParam{8,c},CellParam{9,c},param);
    if c==1
        Par_Array = Par_Array_supp;
    else
        Par_Array = cat(1,Par_Array,Par_Array_supp);
    end
end

%% Randomized
if randomised ==1
    ind = randperm(size(Par_Array,1))';
    if multislice ==1
        ind_minTR = find(Par_Array(:,9)==min(Par_Array(:,9)));
        ind_lowB = find(Par_Array(ind_minTR,7)==min(Par_Array(ind_minTR,7)));
        ind_minTE = find(Par_Array(ind_lowB,8)==min(Par_Array(ind_lowB,8)));
        first_ind = ind_minTE(1,1);
        ind = ind(ind~=first_ind);
        ind = [first_ind; ind];
    else
        if ~isempty(find(Par_Array(:,1)==0))% if short TE experiment
            ind_Short_TE = find(Par_Array(:,1)==0);
            ind_minTE = find(Par_Array(ind_Short_TE,8)==min(Par_Array(ind_Short_TE,8)));
            ind_maxTR = find(Par_Array(ind_minTE,9)==max(Par_Array(ind_minTE,9)));
            first_ind = ind_maxTR(1,1);
            ind = ind(ind~=first_ind);
            ind = [first_ind; ind];
        else                      % find a b0
            ind_B0 = find(Par_Array(:,7)==min(Par_Array(:,7)));
            ind_minTE = find(Par_Array(ind_B0,8)==min(Par_Array(ind_B0,8)));
            ind_maxTR = find(Par_Array(ind_minTE,9)==max(Par_Array(ind_minTE,9)));
            first_ind = ind_maxTR(1,1);
            ind = ind(ind~=first_ind);
            ind = [first_ind; ind];
        end
    end
    Par_Array = Par_Array(ind,:);
else
    ind = 1:size(Par_Array,1);
    if multislice ==1
        ind_minTR = find(Par_Array(:,9)==min(Par_Array(:,9)));
        ind_lowB = find(Par_Array(ind_minTR,7)==min(Par_Array(ind_minTR,7)));
        ind_minTE = find(Par_Array(ind_lowB,8)==min(Par_Array(ind_lowB,8)));
        first_ind = ind_minTE(1,1);
        ind = ind(ind~=first_ind);
        ind = [first_ind; ind];
    end
    Par_Array = Par_Array(ind,:);
end

%% Compute approximation of Frequency
for ind= 1:size(Par_Array,1)
    frequency(1,ind) = calc_freq(Par_Array(ind,3)/1000, Par_Array(ind,4), Par_Array(ind,2));
end

%% Create List
List(:,1) = Par_Array(:,1);% ShortTE
List(:,3) = Par_Array(:,3);% ShapeDuration
Shape_Nlist = Par_Array(:,2) -1; % Shape number between 0 and 5
Shape_Nlist(Shape_Nlist==-1)=0;
List(:,2) = Shape_Nlist;

zeta = acos(sqrt((Par_Array(:,4)*2+1)/3));
g.a = sin(zeta).*Par_Array(:,7);
g.b = sin(zeta).*Par_Array(:,7);
g.c = cos(zeta).*Par_Array(:,7);

alpha = Par_Array(:,5);
beta = Par_Array(:,6);
gamma = zeros(size(List,1));

for n = 1:size(List,1)
    [rotmat,~] = tm_euler_angles2rotmat(alpha(n),beta(n),gamma(n));
    R.xx(n,1) = rotmat(1,1);
    R.xy(n,1) = rotmat(1,2);
    R.xz(n,1) = rotmat(1,3);
    R.yx(n,1) = rotmat(2,1);
    R.yy(n,1) = rotmat(2,2);
    R.yz(n,1) = rotmat(2,3);
    R.zx(n,1) = rotmat(3,1);
    R.zy(n,1) = rotmat(3,2);
    R.zz(n,1) = rotmat(3,3);
end

g.xa = R.xx.*g.a;
g.xb = R.xy.*g.b;
g.xc = R.xz.*g.c;
g.ya = R.yx.*g.a;
g.yb = R.yy.*g.b;
g.yc = R.yz.*g.c;
g.za = R.zx.*g.a;
g.zb = R.zy.*g.b;
g.zc = R.zz.*g.c;
List(:,4) = g.xa;
List(:,5) = g.xb;
List(:,6) = g.xc;
List(:,7) = g.ya;
List(:,8) = g.yb;
List(:,9) = g.yc;
List(:,10) = g.za;
List(:,11) = g.zb;
List(:,12) = g.zc;

List(:,13:14) = Par_Array(:,8:9); % ADD TE TR
List(:,15) = Par_Array(:,7); % Grel

%% Display
Acquisition_time = seconds((sum(sum(List(:,13:14))) +150.*size(List,1))/1000);
Acquisition_time.Format = 'hh:mm';

minTEseq= 20;
figure(1)
set(gcf,'color','w');
sbp(1) = subplot(6,1,1);
plot(List(:,14),'.k')
ylabel('TR (ms)')
ylim([0 ceil(max(List(:,14))/1000)*1000]);
yticks([0 1000 2000 3000 4000 5000])
title(['R1 R2 DOR protocol: ' num2str(size(List,1)) ' experiments, duration = ' char(Acquisition_time) ' hh:mm']);
sbp(2) = subplot(6,1,2);
ShapeTE = 2*List(:,3);
ShapeTE(List(:,1)==0) =0;
plot(List(:,13)+minTEseq+ShapeTE,'.k')
ylabel('TE (ms)')
% yticks([20 30 40 50 60])
sbp(3) = subplot(6,1,3);
EffiFactor_disp = [0 param.EffiFactor];
bvalues_display = 2*param.EffiFactorNorm.*((List(:,15)./(sqrt(param.EffiFactorNorm./(param.EffiFactor(List(:,2)+1))')))*param.Gmax).^2.*(List(:,3)*10^-3).^3.*(gamma1H.^2)*10^-9;
plot(bvalues_display,'.k')
ylabel('b (ms/\mum²)')
yticks([0 1 2 3 4 5])
sbp(4) = subplot(6,1,4);
plot(Par_Array(:,4),'.k')
ylabel('b_{\Delta}')

sbp(5) = subplot(6,1,5);
plot(frequency,'.k')
ylabel('\omega_{cent}/2\pi Hz')

sbp(6) = subplot(6,1,6);
plot(Par_Array(:,5),'.k')
ylim([-pi() pi])
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
ylabel('\Phi')

linkaxes(sbp,'x')

Max_lenght_Shape = max(List(:,3));

%% Write list
% order:IFdiff shape_number,ShapeDuration, Gxa, Gxb, Gxc, Gya, Gyb, Gyc, Gza, Gzb, Gzc, addTE, addTR, grel
list_size = size(List,1);

% Write list
if multislice ==1
    list_file = fopen([out_path, filesep,'DOR_R1R2_',num2str(list_size) '_',num2str(max_shape_order_name) '_ms'],'w');
else
    list_file = fopen([out_path, filesep,'DOR_R1R2_',num2str(list_size) '_',num2str(max_shape_order_name)],'w');
end

fprintf(list_file,'%s\n',['##LIST_SIZE=' num2str(list_size)]);
fprintf(list_file,'%s\n','##LISTDATA');

for n = 1:list_size
    fprintf(list_file,'%i %i %f %f %f %f %f %f %f %f %f %f %f %f %f\n',List(n,1),List(n,2),List(n,3),List(n,4),List(n,5),List(n,6),List(n,7),List(n,8),List(n,9),List(n,10),List(n,11),List(n,12),List(n,13),List(n,14),List(n,15));
end

fprintf(list_file,'%s\n','##END');
fclose(list_file);
disp('list created')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Functions                                          %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phi, theta] = get_dir(pa_path, Ndir)
if Ndir>10
    load(fullfile(pa_path,num2str(Ndir)));
else
    pa_path_filip = split(pa_path,filesep); pa_path_filip = join(pa_path_filip(1:end-1),filesep,1); pa_path_filip = pa_path_filip{1};
    dir_filip =load(fullfile(pa_path_filip,'Elstat_multiopt',['Grad_dirs_multiopt_' num2str(Ndir) '.txt']));
    [phi,theta,~] = cart2sph(dir_filip(:,1),dir_filip(:,2),dir_filip(:,3));
end
end

function [Par_Array] = create_array(pa_path,diff_exp, bval, Ndir,Ndir_spherical,Bdelta,freq_shape,N_b0,TE,TR,param)
gamma1H = 26.75e7;
% bval = [5]; % 10e9 sm-2
% Ndir = 21;
% Ndir_spherical = 4;
% Bdelta = [1 0 -0.5];
% freq_shape = [1 2 3 4 5 6];
% TE = [0];
% TR = [2000 4500];

% Short TE ; Nshape ; Dur Shape ; bDelta ; directions : phi theta ; amplitude ; AddTE ; AddTR
ind_array = 1;
Min_Dur = ceil(((bval'./(2*param.EffiFactorNorm*(param.Gmax^2)*(gamma1H^2)*10^-9)).^(1/3))*10^3);
Min_Dur = max([max(Min_Dur(:));param.Min_Shape_Dur]);
if diff_exp ==1
    Grel = (bval'./(2*param.EffiFactorNorm*(param.Gmax^2).*(Min_Dur*10^-3).^3*(gamma1H^2)*10^-9)).^(1/2);
    %      bval_verif =2*param.EffiFactor.*(Grel.*param.Gmax).^2.*(Min_Dur*10^-3).^3.*(gamma1H.^2)*10^-9;
else
    Grel = 0;
end
if ~isnan(Ndir) && Ndir~=0
    [phi, theta] = get_dir(pa_path, Ndir);
else
    phi =0;  theta =0;
end
if ~isnan(Ndir_spherical) && Ndir_spherical~=0
    [phi_sphere, theta_sphere] = get_dir(pa_path, Ndir_spherical);
else
    phi_sphere =0; theta_sphere =0;
end

if isnan(Ndir_spherical);Ndir_spherical =1;end
if isnan(Ndir);Ndir =1;end

for ind_bval = 1:size(bval,2)
    for ind_Bdelta = 1:size(Bdelta,2)
        if Bdelta(ind_Bdelta)==0
            for ind_dir = 1:Ndir_spherical
                for ind_freq_shape = 1:size(freq_shape,2)
                    if (ind_bval==1 && ind_Bdelta ==1 && ind_dir ==1 && ind_freq_shape ==1)
                        for ind_b0 = 1:N_b0
                            for ind_TE = 1:size(TE,2)
                                for ind_TR = 1:size(TR,2)
                                    Par_Array(ind_array,1) = diff_exp; % Short TE, 1 for no
                                    Par_Array(ind_array,2) = 0; % Nshape
                                    Par_Array(ind_array,3) = Min_Dur; % DurShape
                                    Par_Array(ind_array,4) = 0; % Bdelta
                                    Par_Array(ind_array,5:6) = [0 0]; % directions : phi theta
                                    Par_Array(ind_array,7) = 0; % Amplitude
                                    Par_Array(ind_array,8) = TE(ind_TE); % AddTE
                                    Par_Array(ind_array,9) = TR(ind_TR); % AddTR
                                    ind_array = ind_array +1;
                                end
                            end
                        end
                    end
                    for ind_TE = 1:size(TE,2)
                        for ind_TR = 1:size(TR,2)
                            Par_Array(ind_array,1) = diff_exp; % Short TE, 1 for no
                            Par_Array(ind_array,2) = freq_shape(ind_freq_shape); % Nshape
                            Par_Array(ind_array,3) = Min_Dur; % DurShape
                            Par_Array(ind_array,4) = Bdelta(ind_Bdelta); % Bdelta
                            Par_Array(ind_array,5:6) = [phi_sphere(ind_dir) theta_sphere(ind_dir)]; % directions : phi theta
                            Par_Array(ind_array,7) = Grel(ind_bval).*sqrt(param.EffiFactorNorm./param.EffiFactor(freq_shape(ind_freq_shape)));  % Amplitude
                            Par_Array(ind_array,8) = TE(ind_TE); % AddTE
                            Par_Array(ind_array,9) = TR(ind_TR); % AddTR
                            ind_array = ind_array +1;
                        end
                    end
                end
            end
        else
            for ind_dir = 1:Ndir
                for ind_freq_shape = 1:size(freq_shape,2)
                    if (ind_bval==1 && ind_Bdelta ==1 && ind_dir ==1 && ind_freq_shape ==1)
                        for ind_b0 = 1:N_b0
                            for ind_TE = 1:size(TE,2)
                                for ind_TR = 1:size(TR,2)
                                    Par_Array(ind_array,1) = diff_exp; % Short TE, 1 for no
                                    Par_Array(ind_array,2) = 0; % Nshape
                                    Par_Array(ind_array,3) = Min_Dur; % DurShape
                                    Par_Array(ind_array,4) = 0; % Bdelta
                                    Par_Array(ind_array,5:6) = [0 0]; % directions : phi theta
                                    Par_Array(ind_array,7) = 0; % Amplitude
                                    Par_Array(ind_array,8) = TE(ind_TE); % AddTE
                                    Par_Array(ind_array,9) = TR(ind_TR); % AddTR
                                    ind_array = ind_array +1;
                                end
                            end
                        end
                    end
                    for ind_TE = 1:size(TE,2)
                        for ind_TR = 1:size(TR,2)
                            Par_Array(ind_array,1) = diff_exp; % Short TE, 1 for no
                            Par_Array(ind_array,2) = freq_shape(ind_freq_shape); % Nshape
                            Par_Array(ind_array,3) = Min_Dur; % DurShape  (ind_bval,ind_freq_shape)
                            Par_Array(ind_array,4) = Bdelta(ind_Bdelta); % Bdelta
                            Par_Array(ind_array,5:6) = [phi(ind_dir) theta(ind_dir)]; % directions : phi theta
                            Par_Array(ind_array,7) = Grel(ind_bval).*sqrt(param.EffiFactorNorm./param.EffiFactor(freq_shape(ind_freq_shape))); % Amplitude  (ind_bval,ind_freq_shape)
                            Par_Array(ind_array,8) = TE(ind_TE); % AddTE
                            Par_Array(ind_array,9) = TR(ind_TR); % AddTR
                            ind_array = ind_array +1;
                        end
                    end
                end
            end
        end
    end
end
%% No diff case: diff_exp = 0
if diff_exp==0
    for ind_TE = 1:size(TE,2)
        for ind_TR = 1:size(TR,2)
            Par_Array(ind_array,1) = diff_exp; % Short TE, 1 for no
            Par_Array(ind_array,2) = 0; % Nshape
            Par_Array(ind_array,3) = Min_Dur; % DurShape  (ind_bval,ind_freq_shape)
            Par_Array(ind_array,4) = 0; % Bdelta
            Par_Array(ind_array,5:6) = [0 0]; % directions : phi theta
            Par_Array(ind_array,7) = 0; % Amplitude  (ind_bval,ind_freq_shape)
            Par_Array(ind_array,8) = TE(ind_TE); % AddTE
            Par_Array(ind_array,9) = TR(ind_TR); % AddTR
            ind_array = ind_array +1;
        end
    end
end

Par_Array(isnan(Par_Array))=0;
end

function [frequency] = calc_freq(duration, b_delta, shapeType)
table(:,:,1)=[-1.1986   -1.3667   -0.7574;-2.0613   -2.2706   -1.5604;-3.1023   -3.3595   -2.5102; -4.1938   -4.4730   -3.5705; -5.3047   -5.5965   -4.6669; -6.4251   -6.7248   -5.7793].*1.0e+08;
table(:,:,2)=[0.5568    0.6349    0.3518;0.9576    1.0548    0.7249;1.4412    1.5606    1.1661;1.9482    2.0779    1.6586;2.4643    2.5998    2.1680;2.9848    3.1240    2.6848].*1.0e+07 ;
table(:,:,3)=[-0.8901   -1.0149   -0.5624;-1.5307   -1.6862   -1.1588;-2.3038   -2.4948   -1.8641;-3.1143   -3.3217   -2.6514;-3.9393   -4.1560   -3.4657;-4.7713   -4.9939   -4.2918].*1.0e+05;
table(:,:,4)=[0.5680    0.6476    0.3589;0.9768    1.0760    0.7394;1.4701    1.5920    1.1895;1.9873    2.1196    1.6919;2.5138    2.6520    2.2115;3.0447    3.1867    2.7387].* 1.0e+03;

if b_delta ==0; ind_b_delta = 1;end
if b_delta ==-0.5; ind_b_delta = 2;end
if b_delta ==1; ind_b_delta = 3;end
if b_delta ==0.5; ind_b_delta = 3;end

if shapeType ~=0
    frequency = table(shapeType,ind_b_delta,1)*duration^3 + table(shapeType,ind_b_delta,2)*duration^2 + table(shapeType,ind_b_delta,3)*duration +table(shapeType,ind_b_delta,4);
else
    frequency = 0;
end
end

