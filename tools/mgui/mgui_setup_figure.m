function EG = mgui_setup_figure(EG, c_mode)
% function EG = mgui_setup_figure(EG)
%
% Clear the figure and setup content

if (nargin < 2), c_mode = 2; end

EG.c_mode = c_mode;

% Find figure
if (isfield(EG.handles,'h_fig'))
    figure(EG.handles.h_fig);
else
    error('Must create figure first');
end

if (isfield(EG.handles, 'pix_pos'))
    pix_pos = EG.handles.pix_pos;
end

% Clear the figure and setup content
clf;
set(gcf, ...
    'Color',[1 1 1], ...
    'Unit','Pixels', ...
    'Position', [100 100 1100 710], ...
    'Toolbar', 'none', ...
    'Menubar', 'none', ...
    'Name', ['MDM GUI Version ' EG.conf.version_str], ...
    'CloseRequestFcn', @mgui_close, ...
    'ResizeFcn', @eval_gui_resize, ...
    'NumberTitle','off');

pix_pos.fig = get(gcf,'Position');

% Set up callbacks
set(gcf,'KeyPressFcn', EG.f_callback);
set(gcf,'WindowScrollWheelFcn', EG.f_roi_gui_callback);
set(gcf,'DefaultAxesFontSize', EG.conf.default_font_size)


% Create the panels
ws = 250;
wp = 350;

f_select_pos = @(p)([ 5 5 ws p(4)-10 ]);

h_select = uipanel(...
    'Title','Select data',...
    'FontSize',EG.conf.default_font_size,...
    'BackgroundColor','white',...
    'Units', 'Pixels', ...
    'Position', f_select_pos(pix_pos.fig), ...
    'Tag', EG.t_DATA_PANEL);

pix_pos.data = get(h_select, 'Position');


% -- ROI panel
f_roi_pos = @(p)([ 5+ws+5 5 p(3)-15-ws-wp p(4)-10] );

h_roi = uipanel(...
    'Title', 'Draw VOI',...
    'FontSize',EG.conf.default_font_size,...
    'BackgroundColor','white',...
    'Units', 'pixels', ...
    'Position', f_roi_pos(pix_pos.fig), ...
    'Tag', EG.t_ROI_PANEL);

pix_pos.roi = get(h_roi, 'Position');


% -- ANALYSIS panel
f_analysis_pos = @(p)([ 5+ws+5+p(3)-15-ws-wp 5 wp p(4)-10]); 

h_analysis = uipanel(...
    'Title', 'Analysis',...
    'FontSize',EG.conf.default_font_size,...
    'BackgroundColor','white',...
    'Units', 'pixels', ...
    'Position', f_analysis_pos(pix_pos.fig), ...
    'Tag', EG.t_ANALYSIS_PANEL);

pix_pos.analysis = get(h_analysis, 'Position');



% Local setups
EG = mgui_browse_setup_gui(EG, h_select);
EG = mgui_roi_setup_gui(EG, h_roi);
EG = mgui_analysis_setup_gui(EG, h_analysis);

% Save some handles to the RG structure
EG.handles.pix_pos       = pix_pos;
EG.handles.h_select      = h_select;
EG.handles.h_roi         = h_roi;
EG.handles.h_analysis    = h_analysis;

% Using a struct h_function_list instead of array, since handles are no
% longer doubles
EG.handles.h_function_list.h_roi = h_roi;


% Resizing
    function eval_gui_resize(h, varargin)
        
        if (~exist('EG', 'var')), return; end
        if (~isequal(h, EG.handles.h_fig)), return; end
        if (~isfield(EG, 'handles')), return; end
        if (~isfield(EG.handles, 'h_function_list')), return; end
        
        % Ensure we're working in pixels and get the position
        old_units = get(h,'Units');
        set(h,'Units','pixels');
        p = get(h,'Position');
        
        set(h_select, 'position', f_select_pos(p));
        set(h_roi, 'position', f_roi_pos(p));
        set(h_analysis, 'position', f_analysis_pos(p));
        
        % Restore old units
        set(h,'Units',old_units);
        
    end
end



