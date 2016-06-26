function EG = mgui_analysis_update_panel(EG, do_reload, do_update)
% function EG = mgui_analysis_update_panel(EG, S)

if (nargin < 2), do_reload = 1; end;
if (nargin < 3), do_update = 1; end

EG.roi = msf_ensure_field(EG.roi, 'xps', []);


h_top = EG.handles.h_analysis_top_axes;
h_bottom = EG.handles.h_analysis_bottom_axes;
h_popup  = EG.handles.h_analysis_popup;


% repopulate the model selection dropdown
if (do_reload)
    
    % get current selection
    s = get(h_popup, 'String');
    c = get(h_popup, 'Value');
    try
        current_selection = s{c};
    catch
        current_selection = {};
    end
    
    
    d = dir(fullfile(fileparts(mfilename('fullpath')), '..', '..','models'));
    
    str = {'Overview'};
    for c = 1:numel(d)
        if (d(c).name(1) ~= '.') && (d(c).isdir)
            
            f_name = [d(c).name '_check_xps'];
            
            try
                % currently the check xps functions throw errors
                feval(f_name, EG.roi.xps);
                str{end+1} = d(c).name;

            catch me
                if (0)
                    disp(me.message);
                end
            end
        end
    end
    
    set(h_popup, 'String', str);
    
    value = -1;
    for c = 1:numel(str)
        if (strcmp(str{c}, current_selection))
            value = c;
            break;
        end
    end
    
    if (value > -1)
        set(h_popup, 'value', value);
    end
    
    if (get(h_popup, 'value') > numel(str))
        set(h_popup, 'value', 1);
    end
    
    
    
    value = get(h_popup, 'value');
    switch (value);
        case 1 % show the default overview 1D or 2D
            if (size(EG.roi.I,4) == 1)
                f_plot = @(S, xps) plot_histogram(S, xps, h_top, h_bottom);
            else
                f_plot = @(S, xps) mgui_analysis_plot_overview(S, xps, h_top, h_bottom);
            end
        otherwise
            f_plot = @(S,xps) mgui_analysis_plot(str{value}, S, xps, h_top, h_bottom); 
    end
                
    
    

end

% update the model fits et c
if (do_update)
    
    % pull out signal
    if (msf_isfield(EG, 'roi') && ...
            msf_isfield(EG.roi, 'I') && ...
            msf_isfield(EG.roi, 'I_roi'))
        
        if (~isfield(EG,'analysis') || ~isfield(EG.analysis, 'S'))
            S = zeros(sum(EG.roi.I_roi(:)), size(EG.roi.I,4));
            for c = 1:size(EG.roi.I,4)
                tmp = EG.roi.I(:,:,:,c);
                S(:,c) = tmp(EG.roi.I_roi(:) > 0);
            end
            
            EG.analysis.S = S;
        else
            S = EG.analysis.S;
        end
        
    else
        S = [];
    end
    
    % clear content
    cla(h_top,'reset');
    cla(h_bottom,'reset');
    
    
    xps = EG.roi.xps;
    xps.c_volume = EG.roi.c_volume; % bit nasty, but...
    
    % run analysis and/or plot script
    f_plot(S, xps);
        
end

    function plot_histogram(S, xps, h_top, h_bottom)
        
        hist(h_top, S);
        axis(h_top, 'on');
        
    end


end