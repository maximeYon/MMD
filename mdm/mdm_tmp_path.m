function tmp_path = mdm_tmp_path(do_mkdir)
% function tmp_path = mdm_tmp_path()

% Create a unique temporary directory name to make the code thread safe.
tmp_path = fullfile(tempdir, 'mdm', char(java.util.UUID.randomUUID));

if (do_mkdir)
    mdm_mkdir(tmp_path);
end