DirList = dir;
[DirSize1,DirSize2] = size(DirList);

DirNams = {};
FNams = {};
DirCount = 0;
FileCount = 0;
for ndir = 1:DirSize1
    Namtest = getfield(DirList,{ndir,1},'name');
    if strcmp(Namtest(1),'.') == 0
    if getfield(DirList,{ndir,1},'isdir')
        DirCount = DirCount+1;
        DirNams{DirCount} = Namtest;
    else
        FileCount = FileCount+1;
        FNams{FileCount} = Namtest;
    end
    end
end
