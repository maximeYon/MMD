function [res_fn, tp_fn] = elastix_run_elastix(i_fn, ref_fn, p_fn, o_path)
% function [res_fn, tp_fn] = elastix_run_elastix(i_fn, ref_fn, p_fn, o_path)
%
% Runs elastix. Only for Mac/*nix at the moment. For help, see readme.txt
% in the elastix folder.

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
    error('elastix for windows not implemented');
end

msf_delete({res_fn, tp_fn});
[r, msg] = system(cmd_full);

if (r ~= 0) || (~exist(res_fn, 'file'))
    disp(msg);
    msg = 'If elastix is not installed, check readme.txt in the elastix folder';
    error(sprintf('Could not run ElastiX (%s)\n\n%s', cmd_full, msg));
end

