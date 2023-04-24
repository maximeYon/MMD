function [status, result] = my_FSL_TopUp_Apply(BlipUpImg, BlipDownImg, myFieldMap, nii_fn_out)
% https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/TopupUsersGuide
% 
% This function requires that FSL is installed on your computer. 

% applytopup --imain=my_blipup,my_blipdown --datain=my_parameters --inindex=1,2 --topup=my_field --out=my_good_images

if ispc == 1 % Suppose a WSL installation
    cmd = 'bash -c -i "applytopup --imain=';
    cmd = [cmd '/mnt/' strrep(strrep(lower(BlipUpImg),'\','/'),':','') ',']; % input
    cmd = [cmd '/mnt/' strrep(strrep(lower(BlipDownImg),'\','/'),':','') '']; % input
else
cmd = 'applytopup --imain=';
cmd = [cmd ' "' BlipUpImg '"']; % input
cmd = [cmd ','];
cmd = [cmd ' "' BlipDownImg '"']; % input
end

%% add FSL parameter txt file, use the file created during FieldMap
data_path_pv = split(BlipUpImg,'\');
data_path_pv = join(data_path_pv(1:end-1),filesep,1);
data_path_pv = data_path_pv{1};
% EpiModuleTime=ReadPV360Param([data_path_pv filesep], 'PVM_EpiModuleTime')*1e-3 ;
% PVM_NRepetitions=ReadPV360Param([data_path_pv filesep], 'PVM_NRepetitions');
% PVM_NRepetitions=2;
% FSL_table_pos = repmat([0 1 0 EpiModuleTime],PVM_NRepetitions/2,1);
% FSL_table_neg = repmat([0 -1 0 EpiModuleTime],PVM_NRepetitions/2,1);
% FSL_table = cat(1,FSL_table_pos,FSL_table_neg);
% 
% % Write FSL param to file
% fileID = fopen([data_path_pv filesep 'paramFSL.txt'],'w');
% fprintf(fileID,'%i %i %i %6.4f\n',FSL_table');
% fclose(fileID);

cmd = [cmd ' --datain=' '/mnt/' strrep(strrep([lower(data_path_pv) filesep 'paramFSL.txt'],'\','/'),':','')];

%% add --inindex=1,2
cmd = [cmd ' --inindex=1,2'];

%% add --topup=my_field
cmd = [cmd ' --topup=' '/mnt/' strrep(strrep(lower(myFieldMap),'\','/'),':','')];

%% add output file
% --out=my_good_images

if ispc == 1 % Suppose a WSL installation
    cmd = [cmd ' --out='];
    cmd = [cmd '/mnt/' lower(strrep(strrep(nii_fn_out,'\','/'),':','')) '']; % output 
else
    cmd = [cmd ' --out='];
    cmd = [cmd ' "' nii_fn_out '"']; % output
end

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
% reduced_cmd = cmd{1}(13:end-1);
% clipboard('copy',reduced_cmd)
end

function[my_cnf] =create_cnf(Mean_spatRes, estimated_movement)
my_cnf = ...
{'# (approximate) resolution (in mm) of warp basis for the different sub-sampling levels, default 10';...
['--warpres=' num2str(Mean_spatRes*5) ',' num2str(Mean_spatRes*4) ',' num2str(Mean_spatRes*3.5) ',' num2str(Mean_spatRes*3) ',' num2str(Mean_spatRes*2.5) ',' num2str(Mean_spatRes*1.5) ',' num2str(Mean_spatRes) ',' num2str(Mean_spatRes) ',' num2str(Mean_spatRes)];...
'# sub-sampling scheme, default';...
'--subsamp=1,1,1,1,1,1,1,1,1';...
'# 	FWHM (in mm) of gaussian smoothing kernel, default 8';...
['--fwhm=' num2str(Mean_spatRes*5) ',' num2str(Mean_spatRes*2) ',' num2str(Mean_spatRes*1.5) ',' num2str(Mean_spatRes) ',' num2str(Mean_spatRes*0.75) ',' num2str(Mean_spatRes*0.75) ',' num2str(Mean_spatRes*0.5) ',' num2str(0) ',' num2str(0)];...
'# 	Max # of non-linear iterations, default 5';...
'--miter=5,5,5,5,5,10,10,20,20';...
'# Weight of regularisation, default depending on --ssqlambda and --regmod switches. See user documentation.';...
'--lambda=0.0005,0.0001,0.00001,0.0000015,0.0000005,0.0000005,0.00000005,0.0000000005,0.00000000001';...
'# If set (=1), lambda is weighted by current ssq, default 1';...
'--ssqlambda=1';...
'# Estimate movements if set, default 1 (true)';...
['--estmov=' num2str(estimated_movement) ',' num2str(estimated_movement) ',' num2str(estimated_movement) ',' num2str(estimated_movement) ',' num2str(estimated_movement) ',' num2str(0) ',' num2str(0) ',' num2str(0) ',' num2str(0)];...
'# Minimisation method 0=Levenberg-Marquardt, 1=Scaled Conjugate Gradient, default 0 (LM)';...
'--minmet=0,0,0,0,0,1,1,1,1';...
'# Model for regularisation of warp-field [membrane_energy bending_energy], default bending_energy';...
'--regmod=bending_energy';...
'# Order of spline, 2->Qadratic spline, 3->Cubic spline. Default=3';...
'--splineorder=3';...
'# Precision for representing Hessian, double, or float. Default double';...
'--numprec=double';...
'# Image interpolation model, linear or spline. Default spline';...
'--interp=spline';...
'# 	If set (=1), the images are individually scaled to a common mean, default 0 (false)';...
'--scale=0';...
'# 	If set (=1), the calculations are done in a different grid, default 1 (true)';...
'--regrid=1'};
end
