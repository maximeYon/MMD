function [res_fn, tp_fn] = elastix_run_elastix(i_fn, ref_fn, p_fn, o_path)
% function [res_fn, tp_fn] = elastix_run_elastix(i_fn, ref_fn, p_fn, o_path)
%
% Runs elastix. For help, see readme.txt in the elastix folder.

cmd = 'elastix';
cmd = [cmd ' -f "'   ref_fn  '"']; %#ok<AGROW>
cmd = [cmd ' -m "'   i_fn  '"']; %#ok<AGROW>
cmd = [cmd ' -out "' o_path '"']; %#ok<AGROW>
cmd = [cmd ' -p "'   p_fn  '"']; %#ok<AGROW>

res_fn = fullfile(o_path, 'result.0.nii');
tp_fn  = fullfile(o_path, 'TransformParameters.0.txt');


if (ismac) || (isunix)
    cmd_full = ['/bin/bash --login -c '' ' cmd ' '' '];
else
    cmd_full = cmd;
end

msf_delete({res_fn, tp_fn});
[r, msg] = system(cmd_full);

if (r ~= 0) || (~exist(res_fn, 'file'))
    disp(msg);
    msg = 'If elastix is not installed, check readme.txt in the elastix folder';
    error('%s\n%s\n\n%s\n%s\n%s\n', ...
        'Could not run ElastiX with command:', ...
        cmd_full, ...
        'If elastix is not installed, you will need to install it.', ...
        'Instructions for how to do this is found in the following file', ...
        fullfile(fileparts(mfilename('fullpath')), 'readme.txt'))
end

