function [path, name, ext] = msf_fileparts(fn)

[path, name, ext] = fileparts(fn);

[~,    name, ext2] = fileparts(name);

ext = [ext2 ext];