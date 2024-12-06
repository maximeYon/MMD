    DirList = dir;
    [DirSize1,DirSize2] = size(DirList);
    expno = [];
    for nDir = 1:DirSize1
        expnotest = str2num(getfield(DirList,{nDir,1},'name'));
        if isempty(expnotest) == 0
            expno = [expno expnotest];
        end
    end
    expno = sort(expno,'ascend');
