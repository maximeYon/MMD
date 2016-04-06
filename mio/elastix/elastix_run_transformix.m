function res_fn = elastix_run_transformix(i_fn, t_fn, o_path)
% function res_fn = elastix_run_transformix(i_fn, t_fn, o_path)
%
% Run transformix. Only for Mac/*nix at the moment
% 
% To run, you need to be able to start elastix from your terminal by 
% typing 'transformix'. Do this by adjusting the path in 
% ~/.bashrc according to
%
% export PATH=/usr/local/elastix/bin:/usr/local/elastix/parameters:$PATH
% export DYLD_LIBRARY_PATH=/usr/local/elastix/lib:$DYLD_LIBRARY_PATH

cmd = 'transformix';
cmd = [cmd ' -in "'    i_fn    '"']; 
cmd = [cmd ' -out "'   o_path  '"']; 
cmd = [cmd ' -tp "'    t_fn    '"']; 

res_fn = fullfile(o_path, 'result.nii');

if (ismac) || (isunix)
    cmd_full = ['/bin/bash --login -c '' ' cmd ' '' '];
else
    error('elastix for windows not implemented');
end

msf_delete(res_fn);
[r, msg] = system(cmd_full);

if (r ~= 0) || (~exist(res_fn, 'file'))
    disp(msg);
    error(['Could not run transformix (' cmd_full ')']);
end

    
    
