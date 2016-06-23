function EG = mgui_roi_add_roi(EG, val, bw)
% function EG = mgui_roi_add_roi(EG, val, bw)
%
% This function is closely connected to mgui_roi_format_image
%
% val
%  2 - replace with new ROI
%  1 - add new ROI
% -1 - subtract new ROI

c_dim = EG.roi.c_dim;
c_slice = EG.roi.c_slice;

if (nargin < 3), return; end

I_tmp = mgui_roi_empty(EG);
bw = rot90(bw,-1);
switch (c_dim)
    case 1
        I_tmp(c_slice,:,:) = bw;
    case 2
        I_tmp(:,c_slice,:) = bw;
    case 3
        I_tmp(:,:,c_slice) = bw;
end

if (~isfield(EG.roi, 'I_roi'))
    EG.roi.I_roi = mgui_roi_empty(EG);
    EG.roi.is_updated = 1;
end

% Save for undo clicks
EG.roi.I_roi_undo = EG.roi.I_roi;
set(findobj('Tag', EG.t_ROI_UNDO), 'Enable', 'on');

switch (val)
    case 2 % replace with I_tmp
        EG.roi.I_roi = I_tmp;
        EG.roi.is_updated = 1;
    case 1
        EG.roi.I_roi = EG.roi.I_roi | I_tmp;
        EG.roi.is_updated = 1;
    case -1
        EG.roi.I_roi = EG.roi.I_roi & (~I_tmp);
        EG.roi.is_updated = 1;
    otherwise
        error('strange value');
end

EG.roi = msf_rmfield(EG.roi,'S');

EG = mgui_roi_update_panel(EG, 0, 1);

end