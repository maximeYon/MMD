clear all

% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/';
% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/dicom_100419';
% % % in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/additionalCases';
% out_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/dicom_100419';
% case_dir = dir(fullfile(in_path,'*'));

in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_Spectrum_1/processing_dicom_manmask/acquisition_date_time';
out_path = '/Users/daniel/Dropbox/NMRdata/GE_Spectrum_1/processing_dicom_manmask/acquisition_date_time';
case_dir = dir(fullfile(in_path,'*'));

% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_MSKCC_breast/processing_dicom_manmask/acquisition_date_time';
% out_path = '/Users/daniel/Dropbox/NMRdata/GE_MSKCC_breast/processing_dicom_manmask/acquisition_date_time';
% case_dir = dir(fullfile(in_path,'20*'));

% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_Mahajan_1/processing_dicom_manmask/acquisition_date_time';
% out_path = '/Users/daniel/Dropbox/NMRdata/GE_Mahajan_1/processing_dicom_manmask/acquisition_date_time';
% case_dir = dir(fullfile(in_path,'20*'));
% 
% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/06/acquisition_date_time';
% out_path = '/Users/daniel/Dropbox/NMRdata/Philips_QC/acquisition_date_time';
% case_dir = dir(fullfile(in_path,'*'));

% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/Philips_Sahlgrenska/processing_dicom/acquisition_date_time';
% out_path = '/Users/daniel/Dropbox/NMRdata/Philips_Sahlgrenska/acquisition_date_time';
% case_dir = dir(fullfile(in_path,'*'));

% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_TimSprenger/processing_dicom/2018-12-06_mr450w_index_fix/acquisition_date_time';
% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_TimSprenger/processing_dicom/2018-11-22_Premier_corrected_waveforms/acquisition_date_time';
% out_path = '/Users/daniel/Dropbox/NMRdata/GE_QC/acquisition_date_time';
% case_dir = dir(fullfile(in_path,'*'));

case_names = cell(size(case_dir,1),1);
for ncase = 1:size(case_dir,1)
    case_names{ncase} = case_dir(ncase).name;
end
case_names = case_names(~contains(case_names,'.'));
% case_names = case_names(contains(case_names,'20200303100016'));

msf_mkdir(out_path)
fnam = fullfile(fileparts(out_path),'bootstraps_log.txt');
fid = fopen(fnam, 'w');
fprintf(fid,'%s\n',char(datetime('now'))); 
fprintf(fid,'%s\n',['in_path: ' in_path]); 
fprintf(fid,'%s\n\n',['out_path: ' out_path]); 

for ncase = 1:numel(case_names)
    tic
    case_name = case_names{ncase};
    in_nii_xps_path = fullfile(in_path,case_name,'rwi','nii_2process');
    out_nii_xps_path = fullfile(out_path,case_name,'rwi','nii_xps');
    msf_mkdir(out_nii_xps_path);

    in_fn = fullfile(in_nii_xps_path,'data.nii.gz');
    out_fn = fullfile(out_nii_xps_path,'data.nii.gz');

    if exist(in_fn,'file') == 2 && exist(out_fn,'file') ~= 2
        copyfile(in_fn,out_fn);
        copyfile(mdm_fn_nii2xps(in_fn),mdm_fn_nii2xps(out_fn));
    end

    in_nii_mc_xps_path = fullfile(in_path,case_name,'rwi','res_mc');
    in_fn = fullfile(in_nii_mc_xps_path,'data_mc.nii.gz');
    out_fn = fullfile(out_nii_xps_path,'data_mc.nii.gz');

    if exist(in_fn,'file') == 2 && exist(out_fn,'file') ~= 2
        copyfile(in_fn,out_fn);
        copyfile(mdm_fn_nii2xps(in_fn),mdm_fn_nii2xps(out_fn));
    end
    if exist(mdm_fn_nii2mask(in_fn),'file') == 2 && exist(mdm_fn_nii2mask(out_fn),'file') ~= 2
        copyfile(mdm_fn_nii2mask(in_fn),mdm_fn_nii2mask(out_fn));
    end

    in_bs_path = fullfile(in_path,case_name,'rwi','res_dtd','dtd_res1','bootstraps');
    out_bs_path = fullfile(out_path,case_name,'rwi','dtd','bootstraps');
    Nbs = 0; Nbs_local = 0; Nbs_copied = 0;
    if exist(out_bs_path,'dir') == 7
        bsno_local = msf_getdirno(out_bs_path);
        for nbs = 1:numel(bsno_local)
            out_mfs_fn = fullfile(out_bs_path,num2str(bsno_local(nbs)),'mfs.mat');
            if exist(out_mfs_fn,'file') == 2
                Nbs_local = Nbs_local + 1;
            end
        end
    end
    if exist(in_bs_path,'dir') == 7
        bsno = msf_getdirno(in_bs_path);
        bsno_local = msf_getdirno(out_bs_path);        
        for nbs = 1:numel(bsno)
            in_mfs_fn = fullfile(in_bs_path,num2str(bsno(nbs)),'mfs.mat');
            out_mfs_fn = fullfile(out_bs_path,num2str(bsno(nbs)),'mfs.mat');
            in_ind_fn = fullfile(in_bs_path,num2str(bsno(nbs)),'ind.mat');
            out_ind_fn = fullfile(out_bs_path,num2str(bsno(nbs)),'ind.mat');
            msf_mkdir(fullfile(out_bs_path,num2str(bsno(nbs))));
            if exist(in_mfs_fn,'file') == 2
                Nbs = Nbs + 1;
                if exist(out_mfs_fn,'file') ~= 2
                    copyfile(in_mfs_fn,out_mfs_fn);
                    copyfile(in_ind_fn,out_ind_fn);
                    Nbs_copied = Nbs_copied + 1;
                end
            end
        end
        in_opt_fn = fullfile(in_bs_path,'opt.mat');
        out_opt_fn = fullfile(out_bs_path,'opt.mat');
        if exist(in_opt_fn,'file') == 2 && exist(out_opt_fn,'file') ~= 2
            copyfile(in_opt_fn,out_opt_fn);
        end
    end
    toc_seconds = toc;
    str = [num2str(ncase) ' case: ' case_name '; bootstraps RWI: ' num2str(Nbs) ', local: ' num2str(Nbs_local) ', copied: ' num2str(Nbs_copied) ';  time: ' num2str(toc_seconds,3) ' s'];
    disp(str)
    fprintf(fid,'%s\n',str);
    
end
fclose(fid);

    
