function EG = mgui_create_figure(EG)
% function EG = mgui_create_figure(EG)

% Do not recrete the figure
if (numel(findobj('Tag',EG.t_FIG)) == 0)
    EG.handles.h_fig = figure;
    set(EG.handles.h_fig,'Tag', EG.t_FIG);
end

