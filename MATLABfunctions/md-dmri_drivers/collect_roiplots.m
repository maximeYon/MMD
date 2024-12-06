clear all


nii_paths = cell(0);
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0240306/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0275752/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0333068/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191002085827/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191002090825/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191002184711/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191010133317/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191010133850/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191021111255/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191024170921/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191024171722/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191025191039/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191107094708/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191108140810/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20191111084057/rwi/pmaps';
nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab/20200121173005/rwi/pmaps';

nii_ext = '.nii.gz';
nii_name = 'dtd_s2000';
roi_name = 'roi_lesion';

for nnii = 1:numel(nii_paths)
    nii_path = nii_paths{nnii};
    project_path = fileparts(fileparts(fileparts(nii_path)));
    roiplots_path = fullfile(project_path,'roiplots');
    [~,data_name,~] = fileparts(fileparts(fileparts(nii_path)));

    nii_fn = fullfile(nii_path,[nii_name nii_ext]);
    roi_fn = fullfile(nii_path,[roi_name nii_ext]);

    [I,nii_h] = mdm_nii_read(nii_fn);
    [roi_I,roi_h] = mdm_nii_read(roi_fn);
    roi_I = logical(roi_I);
    roi_projz = squeeze(sum(sum(roi_I,1),2));
    roi_maxz = max(find(roi_projz==max(roi_projz)));

    msf_mkdir(roiplots_path);
    
    in_roiplot_fn = fullfile(fileparts(nii_path),[roi_name '_' nii_name '.pdf']);
    out_roiplot_fn = fullfile(roiplots_path,[data_name '_' roi_name '_' nii_name '.pdf']);
    res = copyfile(in_roiplot_fn,out_roiplot_fn);

    in_histograms_fn = fullfile(fileparts(nii_path),[roi_name '_histograms.pdf']);
    out_histograms_fn = fullfile(roiplots_path,[data_name '_' roi_name '_histograms.pdf']);
    res = copyfile(in_histograms_fn,out_histograms_fn);

    in_sliceplot_fn = fullfile(fileparts(nii_path),'slicecollage',['slice' num2str(roi_maxz) '.pdf']);
    out_sliceplot_fn = fullfile(roiplots_path,[data_name '_slice' num2str(roi_maxz) '.pdf']);
    res = copyfile(in_sliceplot_fn,out_sliceplot_fn);
    
end
