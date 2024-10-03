clear all


project_paths = cell(0);
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_Lund/processing_dicom/acquisition_date_time';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_TimSprenger/processing_dicom/2018-12-06_mr450w_index_fix/acquisition_date_time';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_TimSprenger/processing_dicom/2018-11-22_Premier_corrected_waveforms/acquisition_date_time';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_Spectrum_1/processing_dicom_manmask/acquisition_date_time';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_MSKCC_breast/processing_dicom_manmask/acquisition_date_time';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/GE_Mahajan_1/processing_dicom_manmask/acquisition_date_time';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/dicom_100419';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/additionalCases';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/Szczepankiewicz_DIB_2019_DATA/processing_dicom/acquisition_date_time';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/Philips_Sahlgrenska/processing_dicom/acquisition_date_time';
project_paths{1+numel(project_paths)} = '/Volumes/Matlab/analysis_data_from_FTP_RO/06/acquisition_date_time';

date_str = datestr(datetime('now','TimeZone','local'),'yyyymmdd');
out_fn = ['/Volumes/Matlab/analysis_data_from_FTP_RO/bootstraps_logs/bootstraps_global_' date_str '.xlsx'];

% log_fields = {'AcquisitionDateTime';'PatientName';'PatientID';'ScannerStudyID';'SeriesNumber';'SeriesDescription'};
% log_fields = {'InstitutionName';'Manufacturer';'ManufacturerModelName';'ProductId';'EquipmentUID';
%     'MagneticFieldStrength';'TransmitCoilName';'ReceiveCoil';'ReceiveCoilName';'SoftwareVersion';...
%     'AcquisitionDateTime';'PatientName';'PatientID';'PatientSex';'PatientAge';'BodyPartExamined';...
%     'ScannerStudyID';'SeriesNumber';'SeriesDescription';'SequenceName';'PulseSequenceName';'PulseSequenceDate';'ProtocolName';'ScanOptions';...
%     'EchoTime';'RepetitionTime';'InversionTime';'NumberOfAverages';'NumberOfTemporalPositions';...
%     'PixelSpacing';'ImageDimensionX';'ImageDimensionY';...
%     'SliceThickness';'SpacingBetweenSlices';'LocationsInAcquisition';...
%     'PercentSampling';'PercentPhaseFieldOfView'
%     };
log_fields = {};

log_table = cell(numel(project_paths),1);
for nproject = 1:numel(project_paths)
    project_path = project_paths{nproject};
    case_dir = dir(fullfile(project_path,'*'));
    case_ids = cell(size(case_dir,1),1);
    for ncase = 1:size(case_dir,1)
        case_ids{ncase} = case_dir(ncase).name;
    end
    case_ids = case_ids(~contains(case_ids,'.')); % Remove hidden directories

    emptycell = cell(numel(case_ids),1);
    log_struct.ProjectPath = emptycell;
    log_struct.CaseID = emptycell;
    for nfield = 1:numel(log_fields)
        log_struct.(log_fields{nfield}) = emptycell;
    end
    log_struct.BootstrapsCompleted = emptycell;
%     log_struct.BootstrapsFirstDate = emptycell;
    log_struct.BootstrapsFirstDate = repmat(datetime(0,0,0,0,0,0), numel(case_ids),1);
%     log_struct.BootstrapsLastDate = emptycell;
    log_struct.BootstrapsLastDate = repmat(datetime(0,0,0,0,0,0), numel(case_ids),1);

    for ncase = 1:numel(case_ids)
        case_id = case_ids{ncase};
        log_struct.ProjectPath{ncase} = project_path;
        log_struct.CaseID{ncase} = case_id;
        
        dcmheaders_fn = fullfile(project_path,case_id,'output_from_dicm2nii','dcmHeaders.mat');
        try
            load(dcmheaders_fn);

            NiftiNames = fieldnames(h);    
            for nnifti = 1:numel(NiftiNames)
                NiftiName = NiftiNames{nnifti};
                if numel(NiftiNames)==1, break, end
                if strcmp(h.(NiftiName).AcquisitionDateTime,case_id), break, end
            end
        catch
            warning(['Missing file ' dcmheaders_fn])
        end
                

        for nfield = 1:numel(log_fields)
            try
                log_struct.(log_fields{nfield}){ncase} = h.(NiftiName).(log_fields{nfield});
            catch
                log_struct.(log_fields{nfield}){ncase} = {};
            end
        end
        
        bs_path = fullfile(project_path,case_id,'rwi','res_dtd','dtd_res1','bootstraps');
        bs_dir = dir(fullfile(bs_path,'*','mfs.mat'));        
        log_struct.BootstrapsCompleted{ncase} = numel(bs_dir);
        if numel(bs_dir)>0
            bs_table = struct2table(bs_dir);
            log_struct.BootstrapsFirstDate(ncase) = datetime(min(cellfun(@datetime,bs_table.date)));
            log_struct.BootstrapsLastDate(ncase) = datetime(max(cellfun(@datetime,bs_table.date)));
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

global_log_table_sorted = sortrows(global_log_table,{'BootstrapsLastDate'},{'descend'});

writetable(global_log_table_sorted,out_fn)
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

