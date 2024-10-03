clear all


% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_Spectrum_1/processing_dicom_manmask/acquisition_date_time';
in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_MSKCC_breast/processing_dicom_manmask/acquisition_date_time';
% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_Mahajan_1/processing_dicom_manmask/acquisition_date_time';
% in_path = '/Volumes/Matlab/analysis_data_from_FTP_RO/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/dicom_100419';
out_fn = fullfile(fileparts(in_path),'dcm_info.xlsx');

case_dir = dir(fullfile(in_path,'*'));
case_ids = cell(size(case_dir,1),1);
for ncase = 1:size(case_dir,1)
    case_ids{ncase} = case_dir(ncase).name;
end
case_ids = case_ids(~contains(case_ids,'.')); % Remove hidden directories

% log_fields = {'AcquisitionDateTime';'PatientName';'PatientID';'ScannerStudyID';'SeriesNumber';'SeriesDescription'};
log_fields = {'InstitutionName';'Manufacturer';'ManufacturerModelName';'ProductId';'MagneticFieldStrength';'ReceiveCoilName';'SoftwareVersion';...
    'AcquisitionDateTime';'PatientName';'PatientID';'PatientSex';'PatientAge';'BodyPartExamined';...
    'ScannerStudyID';'SeriesNumber';'SeriesDescription';'PulseSequenceName';'ScanOptions';...
    'EchoTime';'RepetitionTime';'InversionTime';'NumberOfAverages';...
    'SliceThickness';'SpacingBetweenSlices';'LocationsInAcquisition';'NumberOfTemporalPositions';...
    'ImageDimensionX';'ImageDimensionY';'PixelSpacing';'PercentSampling';'PercentPhaseFieldOfView';
    };
emptycell = cell(numel(case_ids),1);

log_struct.CaseID = emptycell;
for nfield = 1:numel(log_fields)
    log_struct.(log_fields{nfield}) = emptycell;
end

for ncase = 1:numel(case_ids)
    case_id = case_ids{ncase};
    log_struct.CaseID{ncase} = case_id;
    dcmheaders_fn = fullfile(in_path,case_id,'output_from_dicm2nii','dcmHeaders.mat');
    load(dcmheaders_fn);

    NiftiNames = fieldnames(h);    
    for nnifti = 1:numel(NiftiNames)
        NiftiName = NiftiNames{nnifti};
        if strcmp(h.(NiftiName).AcquisitionDateTime,case_id) || numel(NiftiNames)==1, break, end
    end
    
    for nfield = 1:numel(log_fields)
        try
            log_struct.(log_fields{nfield}){ncase} = h.(NiftiName).(log_fields{nfield});
        catch
            log_struct.(log_fields{nfield}){ncase} = 'missing';
        end
    end
end

log_table = struct2table(log_struct);
writetable(log_table,out_fn)
%%
% Filename
% Manufacturer
% InstitutionName
% StationName
% StudyDescription
% SeriesDescription
% ManufacturerModelName
% ProductId
% PatientName
% PatientID
% PatientSex
% PatientAge
% BodyPartExamined
% ScanningSequence
% ScanOptions
% SliceThickness
% RepetitionTime
% EchoTime
% InversionTime
% NumberOfAverages
% SpacingBetweenSlices
% EchoTrainLength
% PercentSampling
% PercentPhaseFieldOfView
% PixelBandwidth
% SoftwareVersion
% ProtocolName
% ReceiveCoilName
% InPlanePhaseEncodingDirection
% FlipAngle
% PatientPosition
% PulseSequenceName
% StudyID
% SeriesNumber
% ImagesInAcquisition
% LocationsInAcquisition
% ImageDimensionX
% ImageDimensionY
% PixelSpacing
% AcquisitionDateTime
% NiftiName
% NumberOfTemporalPositions

