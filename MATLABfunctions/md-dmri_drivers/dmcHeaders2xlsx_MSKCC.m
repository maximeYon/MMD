clear all


project_paths = cell(0);
% project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/06/processing_dicom_manmask/acquisition_date_time';
% project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_Spectrum_1/processing_dicom_manmask/acquisition_date_time';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_MSKCC_breast/processing_dicom_manmask/acquisition_date_time';
% project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_Mahajan_1/processing_dicom_manmask/acquisition_date_time';
% project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/dicom_100419';
% project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/additionalCases';

out_fn = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_MSKCC_breast/dcm_info.xlsx';
local_out_fn = '/Users/daniel/Dropbox/NMRdata/GE_MSKCC_breast/dcm_info_mskcc.xlsx';

log_fields = {'AcquisitionDateTime';'PatientName';'PatientID';'PatientAge';'ScannerStudyID';'SeriesNumber';'SeriesDescription';...
     'EchoTime';'RepetitionTime';'InversionTime';'NumberOfTemporalPositions';...
     'PulseSequenceName';'PulseSequenceDate';'ProtocolName';'ScanOptions';...
    'PixelSpacing';'Rows';'Columns';...
    'SliceThickness';'SpacingBetweenSlices';'LocationsInAcquisition';...
   'Pathology';'Location'};
% log_fields = {'InstitutionName';'Manufacturer';'ManufacturerModelName';'ProductId';'EquipmentUID';
%     'MagneticFieldStrength';'TransmitCoilName';'ReceiveCoil';'ReceiveCoilName';'SoftwareVersion';...
%     'AcquisitionDateTime';'PatientName';'PatientID';'PatientSex';'PatientAge';'BodyPartExamined';...
%     'ScannerStudyID';'SeriesNumber';'SeriesDescription';'SequenceName';'PulseSequenceName';'PulseSequenceDate';'ProtocolName';'ScanOptions';...
%     'EchoTime';'RepetitionTime';'InversionTime';'NumberOfAverages';'NumberOfTemporalPositions';...
%     'PixelSpacing';'Rows';'Columns';...
%     'SliceThickness';'SpacingBetweenSlices';'LocationsInAcquisition';...
%     'PercentSampling';'PercentPhaseFieldOfView'
%     };

log_table = cell(numel(project_paths),1);
for nproject = 1:numel(project_paths)
    project_path = project_paths{nproject};
    case_dir = dir(fullfile(project_path,'*'));
    case_ids = cell(size(case_dir,1),1);
    is_hidden = zeros(size(case_dir,1),1);
    for ncase = 1:size(case_dir,1)
        case_ids{ncase} = case_dir(ncase).name;
        is_hidden(ncase) = strcmp(case_ids{ncase}(1),'.');
    end
    case_ids = case_ids(~is_hidden); % Remove hidden directories
%     case_ids = case_ids(~contains(case_ids,'.*')); % Remove hidden directories
%    case_ids = case_ids(~strcmp(case_ids(1),'.')); % Remove hidden directories

    emptycell = cell(numel(case_ids),1);
    log_struct.ProjectPath = emptycell;
    log_struct.CaseID = emptycell;
    for nfield = 1:numel(log_fields)
        log_struct.(log_fields{nfield}) = emptycell;
    end

    for ncase = 1:numel(case_ids)
        case_id = case_ids{ncase};
        log_struct.ProjectPath{ncase} = project_path;
        log_struct.CaseID{ncase} = case_id;
        dcmheaders_fn = fullfile(project_path,case_id,'output_from_dicm2nii','dcmHeaders.mat');
        load(dcmheaders_fn);

        NiftiNames = fieldnames(h);    
        for nnifti = 1:numel(NiftiNames)
            NiftiName = NiftiNames{nnifti};
            if numel(NiftiNames)==1, break, end
            if strcmp(h.(NiftiName).AcquisitionDateTime,case_id), break, end
        end

        for nfield = 1:numel(log_fields)
            try
                log_struct.(log_fields{nfield}){ncase} = h.(NiftiName).(log_fields{nfield});
            catch
                log_struct.(log_fields{nfield}){ncase} = {};
            end
        end
    end
    log_table{nproject} = struct2table(log_struct);
end
%%
global_log_table = log_table{1};
if numel(log_table)>1   
    for ntable = 2:numel(log_table)
        global_log_table = vertcat(global_log_table,log_table{ntable});
    end
end
%%
if exist('out_fn','var')
%     if exist(out_fn,'file')
%         delete(out_fn)
%     end
    writetable(global_log_table,out_fn)
end
%%
if exist('local_out_fn','var')
    if exist(local_out_fn,'file')
        delete(out_fn)
    end
    writetable(global_log_table,local_out_fn)
end

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
% PulseSequenceDate
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

