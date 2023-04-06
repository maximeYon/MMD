% Acquisition protocol for R1 R2 DOR encoding.
%
% Used in Topgaard, Phys. Chem. Chem. Phys. 18, 8545 (2016).
% http://dx.doi.org/10.1039/c5cp07251d

clearvars;
addpath(genpath('/home/maxime/Matlab/MATLABfunctions'));

% Define path for output running and repulsion folders
run_path = fileparts(mfilename('fullpath'));
out_path = '/opt/PV-360.1.1/prog/curdir/maxime/ParaVision/exp/lists/DTD/';
% out_path = 'C:\Users\Administrateur\Mon_Drive\Matlab\MATLABfunctions\md-dmri\acq\bruker\Paravision360\New_Low_level_List';
framework_path = fileparts(fileparts(fileparts(fileparts(run_path))));
pa_path = fullfile(framework_path,'tools','uvec','repulsion_angles');

%% Define parameters
randomised =1; % randomized option

% Gradients and shapes constant
PVM_GradCalConst = 125566;
gamma = 26.75e7;
param.Gmax= PVM_GradCalConst*2*pi/gamma*1e3;
% param.EffiFactor = [0.0053434 0.0026738 0.0016042 0.00083277 0.00049944 0.00031966];
param.EffiFactor = [0.00031966];
% bval = (2*EffiFactor*(Grel.*Gmax)^2*(tau*10^-3)^3)*10^9;
tau = 5;
bval = (2*param.EffiFactor*(1.*param.Gmax)^2*(tau*10^-3)^3)*10^9;
param.Min_Shape_Dur = 5; % minimal duration of shapes in ms


%% 1 Block short TE
diff_exp = 0;
bval = [NaN]; % 10e9 sm-2
Ndir = [NaN];
Ndir_spherical = [NaN];
Bdelta = [NaN];
freq_shape = [NaN];
N_b0 = 0; % Nb_value per TE and TR conditions
TE = [0 10 30];
TR = [4500 100 2000 800];

[Par_Array_supp] = create_array(pa_path,diff_exp, bval, Ndir,Ndir_spherical,Bdelta,freq_shape,N_b0,TE,TR,param);
Par_Array = Par_Array_supp;

%% 2 Block low b-values, spherical (0) only, 4dir, all freq
diff_exp = 1;
bval = [0.69 5.56]; % 10e9 sm-2
Ndir = [15];
Ndir_spherical = 4;
Bdelta = [1 0 -0.5];
freq_shape = [1 2 3 4 5 6];
N_b0 = 3; % Nb_value per TE and TR conditions
TE = [0 20];
TR = [200 800 2000];

[Par_Array_supp] = create_array(pa_path,diff_exp, bval, Ndir,Ndir_spherical,Bdelta,freq_shape,N_b0,TE,TR,param);
Par_Array = cat(1,Par_Array,Par_Array_supp);

%% 3 Block medium b-values, spherical 4 , linear/planar 21, all freq
diff_exp = 1;
bval = [18.76]; % 10e9 sm-2
Ndir = 15;
Ndir_spherical = 4;
Bdelta = [1 0 -0.5];
freq_shape = [1 2 3 4 5 6];
N_b0 = 3; % Nb_value per TE and TR conditions
TE = [0 10];
TR = [2000 4500];

[Par_Array_supp] = create_array(pa_path,diff_exp, bval, Ndir,Ndir_spherical,Bdelta,freq_shape,N_b0,TE,TR,param);
Par_Array = cat(1,Par_Array,Par_Array_supp);

%% 4 Block high b-values, spherical 4 , linear/planar 21, all freq
diff_exp = 1;
bval = [44.49]; % 10e9 sm-2
Ndir = 21;
Ndir_spherical = 4;
Bdelta = [1 0 -0.5];
freq_shape = [1 2 3 4 5 6];
N_b0 = 3; % Nb_value per TE and TR conditions
TE = [0];
TR = [2000 4500];

[Par_Array_supp] = create_array(pa_path,diff_exp, bval, Ndir,Ndir_spherical,Bdelta,freq_shape,N_b0,TE,TR,param);
Par_Array = cat(1,Par_Array,Par_Array_supp);

%% Randomized
if randomised ==1
    ind = randperm(size(Par_Array,1))';
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
    Par_Array = Par_Array(ind,:);
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
sbp(1) = subplot(5,1,1);
plot(List(:,14),'.k')
ylabel('TR (ms)')
ylim([0 ceil(max(List(:,14))/1000)*1000]);
yticks([0 1000 2000 3000 4000 5000])
title(['R1 R2 DOR protocol: ' num2str(size(List,1)) ' experiments, duration = ' char(Acquisition_time) ' hh:mm']);
sbp(2) = subplot(5,1,2);
plot(List(:,13)+minTEseq+2*List(:,3),'.k')
ylabel('TE (ms)')
yticks([20 30 40 50 60])
sbp(3) = subplot(5,1,3);
EffiFactor_disp = [0 param.EffiFactor];
% bval_verif = (2*param.EffiFactor.*(Grel.*param.Gmax).^2.*(Min_Dur*10^-3).^3)*10^9;
% bvalues_display = (2.*EffiFactor_disp(Par_Array(:,2)+1)'.*(List(:,15)*param.Gmax).^2.*(List(:,3)*10^-3).^3)*10^9;
bvalues_display = (2.*EffiFactor_disp(2)'.*(List(:,15)*param.Gmax).^2.*(List(:,3)*10^-3).^3)*10^9;
plot(bvalues_display,'.k')
ylabel('b (ms/\mum²)')
% yticks([0 1 2 3 4 5])
sbp(4) = subplot(5,1,4);
plot(Par_Array(:,4),'.k')
ylabel('b_{\Delta}')
sbp(5) = subplot(5,1,5);
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
list_file = fopen([out_path, filesep,'DOR_R1R2_',num2str(list_size)],'w');

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
% bval = [5]; % 10e9 sm-2
% Ndir = 21;
% Ndir_spherical = 4;
% Bdelta = [1 0 -0.5];
% freq_shape = [1 2 3 4 5 6];
% TE = [0];
% TR = [2000 4500];

% Short TE ; Nshape ; Dur Shape ; bDelta ; directions : phi theta ; amplitude ; AddTE ; AddTR
ind_array = 1;
if diff_exp ==1
    Min_Dur = ceil(((bval'./(2*param.EffiFactor*(param.Gmax^2)*10^9)).^(1/3))./10^-3); % in ms (for 100 % Grad)
    Min_Dur = max(Min_Dur,param.Min_Shape_Dur);
    Grel = (bval'./(2*param.EffiFactor*(param.Gmax^2).*(Min_Dur*10^-3).^3)/10^9).^(1/2);
    %      bval_verif = (2*param.EffiFactor.*(Grel.*param.Gmax).^2.*(Min_Dur*10^-3).^3)*10^9;
else
    Min_Dur = 0;   Grel = 0;
end
if ~isnan(Ndir)
    [phi, theta] = get_dir(pa_path, Ndir);
else
    phi =0;  theta =0;
end
if ~isnan(Ndir_spherical)
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
                                    Par_Array(ind_array,3) = Min_Dur(ind_bval); % DurShape
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
                            Par_Array(ind_array,3) = Min_Dur(ind_bval); % DurShape
                            Par_Array(ind_array,4) = Bdelta(ind_Bdelta); % Bdelta
                            Par_Array(ind_array,5:6) = [phi_sphere(ind_dir) theta_sphere(ind_dir)]; % directions : phi theta
                            Par_Array(ind_array,7) = Grel(ind_bval); % Amplitude
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
                                    Par_Array(ind_array,3) = Min_Dur(ind_bval); % DurShape
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
                            Par_Array(ind_array,3) = Min_Dur(ind_bval); % DurShape  (ind_bval,ind_freq_shape)
                            Par_Array(ind_array,4) = Bdelta(ind_Bdelta); % Bdelta
                            Par_Array(ind_array,5:6) = [phi(ind_dir) theta(ind_dir)]; % directions : phi theta
                            Par_Array(ind_array,7) = Grel(ind_bval); % Amplitude  (ind_bval,ind_freq_shape)
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
Par_Array(isnan(Par_Array))=0;
end
