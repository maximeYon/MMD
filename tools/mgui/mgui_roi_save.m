function EG = mgui_roi_save(EG, roi_filename, data_type)
% function EG = mgui_roi_save(EG, roi_filename, data_type)

if (nargin < 3), data_type = 'uint8'; end

if (~isfield(EG, 'roi')), return; end
if (~isfield(EG.roi,'I_roi')), return; end
if (numel(EG.roi.I_roi(:)) <= 1), return; end

try
    if (~EG.roi.is_updated), return; end
catch
    return;
end

if (nargin < 2)
    if (~isfield(EG.data,'roi'))
        return; % could happen at exits
    else
        roi_filename = mgui_roi_filename(EG);
    end
end

if (numel(roi_filename) == 0), return; end

roi_path = fileparts(roi_filename);

% Create the path
[~,~] = mkdir(roi_path);

% Get header of current document
h = mgui_misc_load_header_as_nifti(EG);

% Flip ROI volume
I_roi = mgui_misc_flip_volume(EG.roi.I_roi, EG.conf.ori, nifti_oricode(h));


if (EG.conf.save_roi_zipped) % try zipping new ROIs
    [roi_path, roi_name, roi_ext] = fileparts(roi_filename);
    roi_fullname = [roi_name roi_ext];
    roi_ext  = roi_fullname((find(roi_fullname == '.', 1, 'first')):end);
    
    switch (lower(roi_ext))
        case '.nii'
            my_delete(roi_filename);
            roi_filename = [roi_filename '.gz'];
    end
    
end

% Save as nifti
switch (data_type)
    case 'uint8'
        nifti_write(uint8(I_roi), roi_filename, h);
    case 'single'
        nifti_write(single(I_roi), roi_filename, h);
end        

end