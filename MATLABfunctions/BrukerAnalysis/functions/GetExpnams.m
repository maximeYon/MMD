DirList = dir;
[DirSize1,DirSize2] = size(DirList);

ExpNam = {};
DirCount = 0;
for ndir = 1:DirSize1
    ExpNamtest = getfield(DirList,{ndir,1},'name');
    if strcmp(ExpNamtest(1),'.') == 0
    if getfield(DirList,{ndir,1},'isdir')
        DirCount = DirCount+1;
        ExpNam{DirCount} = ExpNamtest;
    end
    end
end
