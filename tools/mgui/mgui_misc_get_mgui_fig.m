function h_fig = mgui_misc_get_mgui_fig()
% function h_fig = mgui_misc_get_mgui_fig()

h_fig = get(0,'Children');
ind = cellfun(@(x) strcmp(get(x, 'tag'), 'mgui_FIG'), num2cell(h_fig));
h_fig = h_fig(ind);
