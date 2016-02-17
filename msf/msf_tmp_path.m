function tmp_path = msf_tmp_path(do_mkdir)
% function tmp_path = msf_tmp_path(do_mkdir)

% Create a unique temporary directory name to make the code thread safe.
tmp_path = fullfile(tempdir, 'mdm', char(java.util.UUID.randomUUID));

if (do_mkdir)
    msf_mkdir(tmp_path);
end