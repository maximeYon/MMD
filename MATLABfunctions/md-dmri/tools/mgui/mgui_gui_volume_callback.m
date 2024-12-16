function EG = mgui_gui_volume_callback(EG, ref, hSender, event_info)
% function EG = mgui_gui_volume_callback(EG, ref, hSender)

if (nargin < 4), event_info = []; end

if (isstruct(event_info))
    if (my_isfield(event_info, 'Key'))
        switch (event_info.Key)
            case 'leftarrow'
                EG = mgui_slice_button(EG, ref, -1);
            case 'rightarrow'
                EG = mgui_slice_button(EG, ref, +1);
        end
    end
end

if (ishandle(hSender))
    switch (get(hSender, 'Tag'))
        case ref.tags.SLICE_PLUS
            EG = mgui_slice_button(EG, ref, +1);
        case ref.tags.SLICE_MINUS
            EG = mgui_slice_button(EG, ref, -1);
        case ref.tags.SAG
            EG = mgui_gui_set_dim(EG, ref, 1);
        case ref.tags.COR
            EG = mgui_gui_set_dim(EG, ref, 2);
        case ref.tags.TRA
            EG = mgui_gui_set_dim(EG, ref, 3);
        case ref.tags.SLICE_SLIDER
            EG = mgui_slice_slider(EG, ref);
        case ref.tags.ZOOM_PLUS
            EG = mgui_zoom(EG, ref, +1);
        case ref.tags.ZOOM_MINUS
            EG = mgui_zoom(EG, ref, -1);
    end
end
end
% ----------------------------------------
function EG = mgui_slice_button(EG, r, val)
if (isfield(EG.(r.fn), r.In))    
    if (EG.(r.fn).c_slice + val> size(EG.(r.fn).(r.In), EG.(r.fn).c_dim))
        EG.(r.fn).c_slice = size(EG.(r.fn).(r.In), EG.(r.fn).c_dim);
    elseif (EG.(r.fn).c_slice + val < 1)
        EG.(r.fn).c_slice = 1;
    else
        EG.(r.fn).c_slice = EG.(r.fn).c_slice + val;
    end
    EG = mgui_gui_update_slice_slider(EG, r);
    EG = r.f_callback(EG);      
end
end
% ----------------------------------------
function EG = mgui_slice_slider(EG, r)
EG.(r.fn).c_slice = round(get(findobj('Tag',r.tags.SLICE_SLIDER),'Value'));
EG = r.f_callback(EG);
end
% ----------------------------------------
function EG = mgui_zoom(EG, r, f)
EG = mgui_zoom_format(EG, r);
EG.(r.fn).zoom.factor(EG.(r.fn).c_dim) = EG.(r.fn).zoom.factor(EG.(r.fn).c_dim) * 2^f;
EG = r.f_callback(EG);
end
% ----------------------------------------
function EG = mgui_zoom_format(EG, r)
if (~isfield(EG.(r.fn),'zoom'))
    EG.(r.fn).zoom.factor = ones(1,3);
end
end
