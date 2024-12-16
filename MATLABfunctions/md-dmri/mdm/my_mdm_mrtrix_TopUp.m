function [status, result] = my_mdm_mrtrix_TopUp(nii_fn_in, nii_fn_out, opt)
%
% https://mrtrix.readthedocs.io/en/dev/reference/commands/dwifslpreproc.html
% 
% This function requires that MRTRIX is installed on your computer. Please
% see installation instructions: 
% https://mrtrix.readthedocs.io/en/latest/installation/package_install.html

if ispc == 1 % Suppose a WSL installation
    cmd = 'bash -c -i "dwifslpreproc';
    cmd = [cmd ' /mnt/' lower(nii_fn_in(1)) strrep(nii_fn_in(3:end),'\','/') '']; % input
    cmd = [cmd ' /mnt/' lower(nii_fn_in(1)) strrep(nii_fn_out(3:end),'\','/') '']; % output
else
cmd = 'dwifslpreproc';
cmd = [cmd ' "' nii_fn_in '"']; % input
cmd = [cmd ' "' nii_fn_out '"']; % output
end

%% add options
data_path_pv = split(nii_fn_in,'\');
data_path_pv = join(data_path_pv(1:end-3),filesep,1);
data_path_pv = data_path_pv{1};
EpiModuleTime=ReadPV360Param([data_path_pv filesep], 'PVM_EpiModuleTime')*1e-3 ;
% -rpe_all -pe_dir lr -readout_time 0.66
cmd = [cmd ' -rpe_all -pe_dir 0 -readout_time ' num2str(EpiModuleTime)];

if ~opt.verbose
    cmd = [cmd ' -quiet']; % noise map
end

if opt.do_overwrite
    cmd = [cmd ' -force']; % noise map
end

if ispc == 1 % Suppose a WSL installation
    cmd = [cmd '"'];
    cmd = string(cmd);
end

[status, result] = msf_system(cmd);
