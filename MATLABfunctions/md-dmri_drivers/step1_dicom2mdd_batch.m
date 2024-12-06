% Convert zipped dicom folder to mdd format

clear all

wd = pwd;

% Define paths to dataset folders
datasets_paths = cell(0);
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients_20200130';
datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd19';
datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd18';
datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd17';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd16';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd15';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd14';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd13';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd12';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd11';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd10';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd09';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd08';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd07';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd06';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd05';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd04';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd03';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/Mdd02';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191029night_v1';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Mahajan/MR750w_20191029pm_v1';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/GE/2018-11-22_Premier_corrected_waveforms';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/NIH_Philips';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/MR750_20190521/RWI_volunteer';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191022';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191023';
% datasets_paths{1+numel(datasets_paths)} = '/Users/daniel/Dropbox/NMRdata/Spectrum/MR750w_20191023_pm';

% Prepare options
opt.nii_ext = '.nii.gz';
opt = mdm_opt(opt);

% Loop over datasets
for ndata = 1:numel(datasets_paths)
    datasets_path = datasets_paths{ndata};
    dicoms_path = fullfile(datasets_path,'dicoms');
    niftis_path = fullfile(datasets_path,'niftis');
    mdd_path = fullfile(datasets_path,'mdd');

    % Convert dicoms to niftis
    fnames = mdm_dir2fnames(dicoms_path); 
    dicoms_fn = fullfile(dicoms_path,fnames{1});
    varargout = dicm2nii(dicoms_fn, niftis_path, opt.nii_ext);
    
    load(fullfile(niftis_path,'dcmHeaders.mat'))
    nifti_names = fieldnames(h);
    
    % Loop over niftis
    for nnifti = 1:numel(nifti_names)
        nifti_name = nifti_names{nnifti};
        meta_struct = getfield(h,nifti_name);
        
        % Convert 'epi2extwd' data to mdd format
        if all([strcmp(meta_struct.('PulseSequenceName'),'epi2extwd') ~contains(nifti_name,'ORIG')])
                
            nifti_fn = fullfile(niftis_path,[nifti_name opt.nii_ext]);
            [I,nii_h] = mdm_nii_read(nifti_fn);

            ManufacturerModelName = getfield(getfield(h,nifti_name),'ManufacturerModelName');

            if strcmp(ManufacturerModelName,'DISCOVERY MR750')
                model_str = 'mr750';
            elseif strcmp(ManufacturerModelName,'DISCOVERY MR750w')
                model_str = 'mr450w';
            elseif strcmp(ManufacturerModelName,'SIGNA Premier')
                model_str = 'Premier';
            end

            if contains(nifti_name,'short')
                protocol_str = 'short';
            elseif contains(nifti_name,'intermediate')
                protocol_str = 'intermediate';
            elseif contains(nifti_name,'long')
                protocol_str = 'long';
            end

            xps = mdm_xps_load(['GE_' model_str '_' protocol_str '_xps.mat']);

            % Remove b0 data
            ind = xps.b>1e6;
            I = I(:,:,:,ind);
            xps = mdm_xps_subsample(xps,ind);

            % Save
            nii_name = 'data';
            nii_xps_path = fullfile(mdd_path,nifti_name,'nii_xps');
            msf_mkdir(nii_xps_path);
            nii_fn = fullfile(nii_xps_path, [nii_name opt.nii_ext]);
            xps_fn = mdm_xps_fn_from_nii_fn(nii_fn);
            mdm_nii_write(I, nii_fn, nii_h);
            mdm_xps_save(xps, xps_fn);
        end
    end    
end

