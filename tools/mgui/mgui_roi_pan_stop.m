function mgui_roi_pan_stop(varargin)

% Clear motion functions
h_fig = findobj('Tag','mgui_FIG');
if (~ishandle(h_fig)), return; end
set(h_fig,'WindowButtonMotionFcn', []);
set(h_fig,'WindowButtonUp', []);


% Get the caxis of current axes, then clear senders userdata
hSender = varargin{1};
h = get(hSender, 'CurrentAxes');
s = get(h, 'UserData');
set(h, 'UserData', []);


% Determine the plot active
h_plot = [...
    findobj(get(h, 'Children'), 'Tag', 'mgui_ROI_PROJ_X_PLOT');
    findobj(get(h, 'Children'), 'Tag', 'mgui_ROI_PROJ_Y_PLOT');
    findobj(get(h, 'Children'), 'Tag', 'mgui_ROI_ZOOM_PLOT')];

if (numel(h_plot) ~= 1), return; end
    

% Update GUI: XXX needs better referencing
if (isfield(s,'dx') && isfield(s,'dy'))

    EG = get(h_fig, 'UserData');
    
    switch (get(h_plot, 'Tag'))
        
        case 'mgui_ROI_PROJ_X_PLOT'
            EG.roi.c_slice = EG.roi.c_slice - round(s.dy);
        
        case 'mgui_ROI_PROJ_Y_PLOT'
            EG.roi.c_slice = EG.roi.c_slice + round(s.dx);
    
        case 'mgui_ROI_ZOOM_PLOT'
            EG.roi.zoom.origin{EG.roi.c_dim} = ...
                EG.roi.zoom.origin{EG.roi.c_dim} + [s.dx s.dy];
    end
    
    EG = mgui_roi_update_panel(EG, 0, 1);
    set(h_fig, 'UserData', EG);
end





