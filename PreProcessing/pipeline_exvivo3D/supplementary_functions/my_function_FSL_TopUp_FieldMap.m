function [status, result] = my_function_FSL_TopUp_FieldMap(nii_fn_in, nii_fn_out,Nimages_b0map)
% https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/TopupUsersGuide
% 
% This function requires that FSL is installed on your computer. 

% topup --imain=all_my_b0_images.nii --datain=acquisition_parameters.txt --config=b02b0.cnf --out=my_output

if ispc == 1 % Suppose a WSL installation
    cmd = 'bash -c -i "topup --imain=';
    cmd = [cmd '/mnt/' lower(nii_fn_in(1)) strrep(nii_fn_in(3:end),filesep,'/') '']; % input

else
cmd = 'topup --imain=';
cmd = [cmd nii_fn_in]; % input
end

%% add FSL parameter txt file
data_path_pv = split(nii_fn_in,filesep);
data_path_pv = join(data_path_pv(1:end-3),filesep,1);
data_path_pv = data_path_pv{1};
EpiModuleTime=ReadPV360Param([data_path_pv filesep], 'PVM_EpiModuleTime')*1e-3 ;
PVM_NRepetitions=Nimages_b0map*2;
FSL_table_pos = repmat([0 1 0 EpiModuleTime],PVM_NRepetitions/2,1);
FSL_table_neg = repmat([0 -1 0 EpiModuleTime],PVM_NRepetitions/2,1);
FSL_table = cat(1,FSL_table_pos,FSL_table_neg);

% Write FSL param to file
fileID = fopen([data_path_pv filesep 'paramFSL.txt'],'w');
fprintf(fileID,'%i %i %i %6.4f\n',FSL_table');
fclose(fileID);

if ispc == 1 % Suppose a WSL installation
cmd = [cmd ' --datain=' '/mnt/' strrep(strrep([lower(data_path_pv) filesep 'paramFSL.txt'],filesep,'/'),':','')];
else
cmd = [cmd  ' --datain=' data_path_pv filesep 'paramFSL.txt'];
end

% %% add FSL parameter cnf file
% SpatDim = ReadPV360Param([data_path_pv filesep], 'PVM_SpatDimEnum');
% PVM_SpatResol=ReadPV360Param([data_path_pv filesep], 'PVM_SpatResol');
% PVM_SpatResol_phase = PVM_SpatResol(2:end);
% estimated_movement = 0; %  default value: 100
[my_cnf] =create_cnf();
% Write FSL cnf to file
fileID = fopen([data_path_pv filesep 'confFSL.cnf'],'w');
for indL = 1:size(my_cnf,1)
fprintf(fileID,'%s\n',my_cnf{indL,1});
end
fclose(fileID);

if ispc == 1 % Suppose a WSL installation
cmd = [cmd ' --config=' '/mnt/' strrep(strrep([lower(data_path_pv) filesep 'confFSL.cnf'],filesep,'/'),':','')];
else
cmd = [cmd  ' --config=' data_path_pv filesep 'confFSL.cnf'];
end
%% add output file
% --out=my_topup_results --fout=my_field --iout=my_unwarped_images
nii_fn_out_path = split(nii_fn_out,filesep);
nii_fn_out_path = join(nii_fn_out_path(1:end-1),filesep,1);
nii_fn_out_path = [nii_fn_out_path{1} filesep 'fsl'];

if ispc == 1 % Suppose a WSL installation
    cmd = [cmd ' --out='];
    cmd = [cmd '/mnt/' lower(strrep(strrep(nii_fn_out_path,filesep,'/'),':','')) '']; % output 
    cmd = [cmd ' --fout='];
    cmd = [cmd '/mnt/' lower(strrep(strrep(nii_fn_out_path,filesep,'/'),':','')) '']; % output 
    cmd = [cmd ' --iout='];
    cmd = [cmd '/mnt/' lower(strrep(strrep(nii_fn_out,filesep,'/'),':','')) '']; % output 
else
    cmd = [cmd ' --out='];
    cmd = [cmd nii_fn_out_path]; % output
    cmd = [cmd ' --fout='];
    cmd = [cmd nii_fn_out_path]; % output
    cmd = [cmd ' --iout='];
    cmd = [cmd nii_fn_out]; % output
end

% %% test verbose
% cmd = [cmd ' --verbose']; % output

%% test paralelized
cmd = [cmd ' --nthr=''24''']; 

if ispc == 1 % Suppose a WSL installation
    cmd = [cmd '"'];
    cmd = string(cmd);
end

[status, result] = msf_system(cmd);
%% If that doesn't work check that the ~/.bashrc file contains
% FSLDIR=/usr/local/fsl
% PATH=${FSLDIR}/bin:${PATH}
% export FSLDIR PATH
% . ${FSLDIR}/etc/fslconf/fsl.sh


%% Alternative, copy to clipboard
reduced_cmd = cmd{1}(13:end-1);
clipboard('copy',reduced_cmd)
end

function[my_cnf] =create_cnf()
my_cnf = ...
{'# (approximate) resolution (in mm) of warp basis for the different sub-sampling levels, default 10';...
'--warpres=20,16,14,12,10,6,4,4,4,2';...
% '--warpres=20,16,14,12,10,6,4,4,4';...
'# sub-sampling scheme, default';...
'--subsamp=1,1,1,1,1,1,1,1,1,1';...
% '--subsamp=1,1,1,1,1,1,1,1,1';...
'# 	FWHM (in mm) of gaussian smoothing kernel, default 8';...
'--fwhm=8,6,4,3,3,2,1,0,0,0';...
% '--fwhm=8,6,4,3,3,2,1,0,0';...
'# 	Max # of non-linear iterations, default 5';...
% '--miter=5,5,5,5,5,8,8,15,15';...
'--miter=5,5,5,5,5,8,8,8,8,8';...
'# Weight of regularisation, default depending on --ssqlambda and --regmod switches. See user documentation.';...
'--lambda=0.005,0.001,0.0001,0.000015,0.00000005,0.000000005,0.0000000005,0.000000000005,0.0000000000001,0.0000000000001';...
% '--lambda=0.005,0.001,0.0001,0.000015,0.00000005,0.000000005,0.0000000005,0.000000000005,0.0000000000001';...
'# If set (=1), lambda is weighted by current ssq, default 1';...
'--ssqlambda=1';...
'# Model for regularisation of warp-field [membrane_energy bending_energy], default bending_energy';...
'--regmod=bending_energy';...
'# Estimate movements if set, default 1 (true)';...
'--estmov=0,0,0,0,0,0,0,0,0,0';...
% '--estmov=0,0,0,0,0,0,0,0,0';...
'# Minimisation method 0=Levenberg-Marquardt, 1=Scaled Conjugate Gradient, default 0 (LM)';...
'--minmet=0,0,0,0,0,1,1,1,1,1';...
% '--minmet=0,0,0,0,0,1,1,1,1';...
'# Order of spline, 2->Qadratic spline, 3->Cubic spline. Default=3';...
'--splineorder=3';...
'# Precision for representing Hessian, double, or float. Default double';...
'--numprec=double';...
'# Image interpolation model, linear or spline. Default spline';...
'--interp=spline';...
'# 	If set (=1), the images are individually scaled to a common mean, default 0 (false)';...
'--scale=1'};...
end

% %% Orig
% function[my_cnf] =create_cnf()
% my_cnf = ...
% {'# (approximate) resolution (in mm) of warp basis for the different sub-sampling levels, default 10';...
% % '--warpres=20,16,14,12,10,6,4,4,4,2';...
% '--warpres=20,16,14,12,10,6,4,4,4';...
% '# sub-sampling scheme, default';...
% % '--subsamp=1,1,1,1,1,1,1,1,1,1';...
% '--subsamp=1,1,1,1,1,1,1,1,1';...
% '# 	FWHM (in mm) of gaussian smoothing kernel, default 8';...
% % '--fwhm=8,6,4,3,3,2,1,0,0,0';...
% '--fwhm=8,6,4,3,3,2,1,0,0';...
% '# 	Max # of non-linear iterations, default 5';...
% '--miter=5,5,5,5,5,8,8,15,15';...
% % '--miter=5,5,5,5,5,8,8,8,8,8';...
% '# Weight of regularisation, default depending on --ssqlambda and --regmod switches. See user documentation.';...
% % '--lambda=0.005,0.001,0.0001,0.000015,0.00000005,0.000000005,0.0000000005,0.000000000005,0.0000000000001,0.0000000000001';...
% '--lambda=0.005,0.001,0.0001,0.000015,0.00000005,0.000000005,0.0000000005,0.000000000005,0.0000000000001';...
% '# If set (=1), lambda is weighted by current ssq, default 1';...
% '--ssqlambda=1';...
% '# Model for regularisation of warp-field [membrane_energy bending_energy], default bending_energy';...
% '--regmod=bending_energy';...
% '# Estimate movements if set, default 1 (true)';...
% % '--estmov=0,0,0,0,0,0,0,0,0,0';...
% '--estmov=0,0,0,0,0,0,0,0,0';...
% '# Minimisation method 0=Levenberg-Marquardt, 1=Scaled Conjugate Gradient, default 0 (LM)';...
% % '--minmet=0,0,0,0,0,1,1,1,1,1';...
% '--minmet=0,0,0,0,0,1,1,1,1';...
% '# Order of spline, 2->Qadratic spline, 3->Cubic spline. Default=3';...
% '--splineorder=3';...
% '# Precision for representing Hessian, double, or float. Default double';...
% '--numprec=double';...
% '# Image interpolation model, linear or spline. Default spline';...
% '--interp=spline';...
% '# 	If set (=1), the images are individually scaled to a common mean, default 0 (false)';...
% '--scale=1'};...
% end
