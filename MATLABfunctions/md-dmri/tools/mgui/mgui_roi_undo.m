function EG = mgui_roi_undo(EG)
% function EG = mgui_roi_undo(EG)

% Save for undo clicks
EG.roi.I_roi = EG.roi.I_roi_undo;
EG.roi.I_roi_undo = [];
set(findobj('Tag', EG.t_ROI_UNDO), 'Enable', 'off');

EG.roi = msf_rmfield(EG.roi,'S');
EG = mgui_roi_update_panel(EG, 0, 1);

