function I_roi = mgui_roi_empty(EG)
% function I_roi = mgui_roi_empty(EG)

if (~isfield(EG,'roi') || (~isfield(EG.roi,'I')))
    I_roi = [];
elseif (ndims(EG.roi.I) <= 3)
    I_roi = zeros(size(EG.roi.I));
elseif (ndims(EG.roi.I) == 4) && (size(EG.roi.I,1) == 3)
    a = size(EG.roi.I);
    I_roi = zeros(a(2:4));
elseif (ndims(EG.roi.I) == 4)
    a = size(EG.roi.I);
    I_roi = zeros(a(1:3));
else
    error('stop');
end
