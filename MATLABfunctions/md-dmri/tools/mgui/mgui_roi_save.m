function EG = mgui_roi_save(EG, roi_filename, data_type)
% function EG = mgui_roi_save(EG, roi_filename, data_type)

if (nargin < 3), data_type = 'uint8'; end

if (~isfield(EG, 'roi')), return; end
if (~isfield(EG.roi,'I_roi')), return; end
if (numel(EG.roi.I_roi(:)) <= 1), return; end
if (numel(roi_filename) == 0), return; end

try
    if (~EG.roi.is_updated), return; end
catch
    return;
end


% Create the path
roi_path = msf_fileparts(roi_filename);
[~,~] = mkdir(roi_path);

% Get header of current document
h = EG.roi.header;

% Flip ROI volume
I_roi = mgui_misc_flip_volume(EG.roi.I_roi, EG.conf.ori, mdm_nii_oricode(h));

% Save as nifti
switch (data_type)
    case 'uint8'
        mdm_nii_write(uint8(I_roi), roi_filename, h);
    case 'single'
        mdm_nii_write(single(I_roi), roi_filename, h);
end        

