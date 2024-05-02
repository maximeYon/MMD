function [hdr,ext] = mgui_contrast_load_header(EG, select)
% function hdr = mgui_contrast_load_header(EG, select)

% XXX: There's a potential problem with the flipping, which is done only in
% nifti_read but not in here. Code should probably be better structures in
% this aspect /MN, 2016

if (nargin < 2), select = EG.select; end

% Get filename
[filename, ext] = mgui_contrast_filename(EG, select);

% Search for .gz files if the orignal file is not present
if (~exist(filename, 'file') && exist([filename '.gz'], 'file'))
    filename = [filename '.gz'];
end


% Load the image volume, if the file exists
hdr = [];
if (exist(filename, 'file'))
    switch (lower(ext))
        case '.nii'
            % make sure to duplicate parts of this in mgui_misc_load_volume,
            % in order to avoid reading the header twice
            hdr = mdm_nii_read_header(filename);
            hdr.my_hdr.ori = mdm_nii_oricode(hdr);
           
                        
        otherwise
            error('not supported');
    end
end


