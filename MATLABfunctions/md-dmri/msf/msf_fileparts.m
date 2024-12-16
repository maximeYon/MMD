function [path, name, ext] = msf_fileparts(fn)
% function [path, name, ext] = msf_fileparts(fn)
%
% Just as the built-in 'fileparts', except that this run fileparts twice in
% order to count '.nii.gz' as one extension

[path, name, ext] = fileparts(fn);

[~,    name, ext2] = fileparts(name);

ext = [ext2 ext];