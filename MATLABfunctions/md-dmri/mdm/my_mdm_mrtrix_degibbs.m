function [status, result] = my_mdm_mrtrix_degibbs(nii_fn_in, nii_fn_out, opt)
% function [status, result] = mdm_mrtrix_degibbs(nii_fn_in, nii_fn_out, opt)
%
% https://mrtrix.readthedocs.io/en/latest/reference/commands/mrdegibbs.html
% 
% This function requires that MRTRIX is installed on your computer. Please
% see installation instructions: 
% https://mrtrix.readthedocs.io/en/latest/installation/package_install.html

if nargin < 3
    opt = mdm_mrtrix_opt();
end

if ispc == 1 % Suppose a WSL installation
    cmd = 'bash -c -i "mrdegibbs';
    cmd = [cmd ' /mnt/' lower(nii_fn_in(1)) strrep(nii_fn_in(3:end),'\','/') '']; % input
    cmd = [cmd ' /mnt/' lower(nii_fn_in(1)) strrep(nii_fn_out(3:end),'\','/') '']; % output
else
cmd = 'mrdegibbs';
cmd = [cmd ' "' nii_fn_in '"']; % input
cmd = [cmd ' "' nii_fn_out '"']; % output
end


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
