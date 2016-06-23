function EG = mgui_roi_set_function(EG, c_function)
% function EG = mgui_roi_set_function(EG, c_function)

EG.roi.function = c_function;

tags = {EG.t_ROI_PLUS, EG.t_ROI_MINUS};


for c = 1:numel(tags)
    
    if (c == c_function)
        bg_color = [36 40 215] / 256;
        color = [1 1 1];
        fw = 'bold';
    else
        bg_color = [0 0 0] + 0.94;
        color = [0 0 0];
        fw = 'normal';
    end
        
    
    set(EG.tags.(tags{c}), ...
        'BackgroundColor', bg_color, ...
        'ForegroundColor', color, ...
        'FontWeigh', fw);
    
end

