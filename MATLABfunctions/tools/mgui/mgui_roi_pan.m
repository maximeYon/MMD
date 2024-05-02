function mgui_roi_pan(varargin)
% function mgui_roi_pan(varargin)


% Get userdata and current point from sender
hSender = varargin{1};
h = get(hSender, 'CurrentAxes');
p = get(hSender, 'CurrentPoint');
s = get(h, 'UserData');
c = get(h, 'Children');

% Determine which plot that is active
h_plot = [];
for c_child = 1:numel(c)
    switch(get(c(c_child), 'Tag'))
        case {'mgui_ROI_PROJ_X_PLOT', 'mgui_ROI_PROJ_Y_PLOT', 'mgui_ROI_ZOOM_PLOT'}
            h_plot = c(c_child); 
            break;
    end
end

% Ensure that only one plot is present
if (numel(h_plot) ~= 1), return; end

% Extract useful information
if (~isfield(s,'p0')),      s.p0 = p; end
if (~isfield(s,'xdata0')),  s.xdata0 = get(h_plot, 'XData');  end
if (~isfield(s,'ydata0')),  s.ydata0 = get(h_plot, 'YData');  end
if (~isfield(s,'ps')),      s.ps = get(h, 'position'); s.ps = s.ps(3:4); end
if (~isfield(s,'ax')),      s.ax = axis(h); s.ax = s.ax([2 4]); end
if (~isfield(s,'as')),      s.as = get(h, 'plotboxaspectratio'); end

s.dx = +(p(1) - s.p0(1)) / s.ps(1) * s.ax(1);
s.dy = -(p(2) - s.p0(2)) / s.ps(2) * s.ax(2)  * (s.as(1)/s.as(2));

% Handle each object separately
switch (get(h_plot, 'Tag'))
    case 'mgui_ROI_PROJ_X_PLOT'
        set(h_plot, 'YData', s.ydata0 + s.dy);
    case 'mgui_ROI_PROJ_Y_PLOT'
        set(h_plot, 'XData', s.xdata0 + s.dx);
    case 'mgui_ROI_ZOOM_PLOT'
        set(h_plot, 'XData', s.xdata0 + s.dx);
        set(h_plot, 'YData', s.ydata0 + s.dy);
end

% Store the settings
set(h, 'UserData', s);

