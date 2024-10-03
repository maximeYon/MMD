clear all

in_path = '/Users/daniel/Dropbox/NMRdata/MSKCC/RWI_Matlab';
out_path = '/Users/daniel/Dropbox/NMRdata/GE_MSKCC_breast/processing_dicom_manmask/rois_DT';
case_dir = dir(fullfile(in_path,'20*'));

case_names = cell(size(case_dir,1),1);
for ncase = 1:size(case_dir,1)
    case_names{ncase} = case_dir(ncase).name;
end

msf_mkdir(out_path)
in_roi_name = 'roi_lesion.nii.gz';
out_roi_name = 'lesion.nii.gz';


for ncase = 1:numel(case_names)
    case_name = case_names{ncase};
    in_fn = fullfile(in_path,case_name,'rwi','pmaps',in_roi_name);
    out_fn = fullfile(out_path,case_name,out_roi_name);
    msf_mkdir(fileparts(out_fn));

    if exist(in_fn,'file') == 2 && exist(out_fn,'file') ~= 2
        copyfile(in_fn,out_fn);
        disp(case_name)
    end
    
end

    
