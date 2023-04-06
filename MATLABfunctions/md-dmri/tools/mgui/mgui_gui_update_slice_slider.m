function EG = mgui_gui_update_slice_slider(EG, r)
% function EG = mgui_gui_update_slice_slider(EG, r)

if (isfield(EG.(r.fn),r.In))
    set(EG.tags.(r.tags.SLICE_SLIDER), ...
        'Min', 0.9, ...
        'Max', size(EG.(r.fn).(r.In), EG.(r.fn).c_dim), ...
        'Value', EG.(r.fn).c_slice);
end

