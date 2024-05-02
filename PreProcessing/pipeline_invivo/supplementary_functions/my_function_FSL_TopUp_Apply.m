function [status, result] = my_function_FSL_TopUp_Apply(BlipUpImg, BlipDownImg, myFieldMap, nii_fn_out)
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
cmd = [cmd BlipUpImg]; % input
cmd = [cmd ','];
cmd = [cmd BlipDownImg]; % input
end

%% add FSL parameter txt file, use the file created during FieldMap
data_path_pv = split(BlipUpImg,filesep);
data_path_pv = join(data_path_pv(1:end-3),filesep,1);
data_path_pv = data_path_pv{1};
EpiModuleTime=ReadPV360Param([data_path_pv filesep], 'PVM_EpiModuleTime')*1e-3 ;
PVM_NRepetitions=2;
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
%% add --inindex=1,2
cmd = [cmd ' --inindex=1,2'];

%% add --topup=my_field
if ispc == 1 % Suppose a WSL installation
cmd = [cmd ' --topup=' '/mnt/' strrep(strrep(lower(myFieldMap),'\','/'),':','')];
else
cmd = [cmd ' --topup=' myFieldMap];
end
%% add output file
% --out=my_good_images

if ispc == 1 % Suppose a WSL installation
    cmd = [cmd ' --out='];
    cmd = [cmd '/mnt/' lower(strrep(strrep(nii_fn_out,'\','/'),':','')) '']; % output 
else
    cmd = [cmd ' --out='];
    cmd = [cmd nii_fn_out]; % output
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
