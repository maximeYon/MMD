function dirno = msf_getdirno(dir_path)
% function dirno = msf_getdirno(dir_path)

dirlist = dir(dir_path);
[Ndir,~] = size(dirlist);
dirno = [];
for ndir = 1:Ndir
    dirnotest = str2num(getfield(dirlist,{ndir,1},'name'));
    if ~isempty(dirnotest)
        dirno = [dirno dirnotest];
    end
end
dirno = sort(dirno,'ascend');

