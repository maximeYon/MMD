function EG = mgui_waitbar(EG, value, message, s_panel)
% function EG = mgui_waitbar(EG, value, message)
%
% value = 0, create new waitbar
% value > 0, update waitbar, if present
% value > 1, close waitbar
%
% Intended use: Create waitbar outside updating function, if waitbar is
% desired.

if (EG.conf.disable_waitbar)
    return;
end

if (nargin < 2), return; end
if (nargin < 3), message = []; end
if (nargin < 4), s_panel = []; end

if (isempty(message)), message = ''; end
if (isempty(s_panel)), s_panel = ''; end
        
if (value < 0)
    
    error('Value must be above zero');
    
elseif (value == 0) % Create new waitbar
    
    % Close current waitbar, if there is one present
    if (isfield(EG,'handles') && isfield(EG.handles,'h_waitbar'))
        EG = mgui_waitbar(EG, inf);
    end
    
    % Store timing of call to waitbar; only open waitbar if there are long
    % opening times
    if (numel(s_panel) > 0)
        EG.waitbar.(s_panel).tic = tic;
        
        if (~isfield(EG.waitbar.(s_panel), 'toc'))
            EG.waitbar.(s_panel).toc = [];
        end
        
        n = 5;
        
        if (numel(EG.waitbar.(s_panel).toc) < n)
            do_open = 1;
        elseif (mean(EG.waitbar.(s_panel).toc( (end-n+1):end )) > 0.5)
            do_open = 1;
        else
            do_open = 0;
        end
        
    else
        do_open = 1;
    end
    
    % Start a new waitbar
    if (do_open)
        EG.handles.h_waitbar = waitbar(value, message, ...
            'WindowStyle', 'Modal');
        drawnow; pause(0.05); % prevents Matlab from hanging?
    end
    
    
elseif (value < 1) % Update waitbar if there already is a handle
    
    if (isfield(EG, 'handles') && isfield(EG.handles,'h_waitbar') ...
            && ishandle(EG.handles.h_waitbar))
        waitbar(value, EG.handles.h_waitbar, message);
        drawnow; pause(0.05); % prevents Matlab from hanging?
    end
    
else % Close waitbar if value is equal to or above unity
    
    if (isfield(EG.handles,'h_waitbar') && ishandle(EG.handles.h_waitbar))
        close(EG.handles.h_waitbar);
        EG.handles = rmfield(EG.handles,'h_waitbar');
    end
    
    if (numel(s_panel) > 0)
        EG.waitbar.(s_panel).toc(end+1) = toc(EG.waitbar.(s_panel).tic);
    end
    
end
