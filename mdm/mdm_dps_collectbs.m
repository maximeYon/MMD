function bs_dps = mdm_dps_collectbs(method, bs_path, opt)
% function bs_dps = mdm_dps_collectbs(method, bs_path, opt)
%

bsno = msf_getdirno(bs_path);        
bs_dps = cell(numel(bsno),1);
parfor nbs = 1:numel(bsno)
    mfs_fn   = fullfile(bs_path,num2str(bsno(nbs)),'mfs.mat');
    if exist(mfs_fn,'file')==2

        mfs = mdm_mfs_load(mfs_fn);
        m = mfs.m;

        if strcmp(method,'dtr2d')
            bs_dps{nbs} = dtr2d_4d_fit2param(m, [], opt);
        elseif strcmp(method,'dtr1d')
            bs_dps{nbs} = dtr1d_4d_fit2param(m, [], opt);
        elseif strcmp(method,'dtr1r2d')
            bs_dps{nbs} = dtr1r2d_4d_fit2param(m, [], opt);
        elseif strcmp(method,'dtd')
            bs_dps{nbs} = dtd_4d_fit2param(m, [], opt);
        end
        
        bs_dps{nbs}.nii_h = mfs.nii_h;
                
    end
end

