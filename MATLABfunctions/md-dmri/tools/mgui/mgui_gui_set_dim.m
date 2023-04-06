function EG = mgui_gui_set_dim(EG, r, c_dim)
% function EG = mgui_gui_set_dim(EG, r, c_dim)

% Save original dimension
if (isfield(EG.(r.fn),'c_dim'))
    c_dim_o = EG.(r.fn).c_dim; 
else
    c_dim_o = 1; 
end
EG.(r.fn).c_dim = c_dim;

if (isfield(EG.(r.fn), r.In))
    
    % Remember slices between dims
    if (~isfield(EG.(r.fn),'c_slice_save'))
        for c = 1:3 % init by looping over dimensions
            EG.(r.fn).c_slice_save(c) = round(size(EG.(r.fn).(r.In), c) / 2);
        end
    else
        EG.(r.fn).c_slice_save(c_dim_o) = EG.(r.fn).c_slice;
    end

    % Check number of slices
    if (ndims(EG.(r.fn).(r.In)) <= 3)
        EG.(r.fn).n_slice = size(EG.(r.fn).(r.In), EG.(r.fn).c_dim);
    elseif (ndims(EG.(r.fn).(r.In)) == 4) && (size(EG.(r.fn).(r.In),1) == 3)
        EG.(r.fn).n_slice = size(EG.(r.fn).(r.In), EG.(r.fn).c_dim + 1);
    elseif (ndims(EG.(r.fn).(r.In)) == 4)
        EG.(r.fn).n_slice = size(EG.(r.fn).(r.In), EG.(r.fn).c_dim);
    else
        error('strange image');
    end
    
    % Set slice and update GUI
    if (numel(EG.(r.fn).c_slice_save) == 2) && (EG.(r.fn).c_dim == 3)
        EG.(r.fn).c_slice = 1;
    else
        EG.(r.fn).c_slice = EG.(r.fn).c_slice_save(EG.(r.fn).c_dim);
    end
    EG = mgui_gui_update_slice_slider(EG, r);
    EG = r.f_callback(EG);
else
    % Clear image
    EG = r.f_callback(EG);
end

