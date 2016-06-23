function eval_gui_close(varargin)
% function eval_gui_close(varargin)

EG = mgui_define_tags;

if (ishandle(findobj('Tag',EG.t_FIG)))
    EG = get(findobj('Tag',EG.t_FIG), 'UserData');    
    mgui_roi_save(EG);
    clear EG;
end

delete(gcf);

return;
