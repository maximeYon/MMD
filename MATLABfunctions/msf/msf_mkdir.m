function msf_mkdir(folder_path)
% function msf_mkdir(folder_path)
%
% Recursively makes folder_path

if (isempty(folder_path)), return; end

if (isunix)
    
    try
        %eval(['!mkdir -p ' folder_path]);
        eval(['!mkdir -p $"' folder_path '"']); %Fix to work with spaces in path. DT 20200201 
    catch 
        disp('could not make directory');
    end
    
else
    
    [status, ~, ~] = mkdir(folder_path);
    
    if (~status)
        error('Could not make directory');
    end
    
end