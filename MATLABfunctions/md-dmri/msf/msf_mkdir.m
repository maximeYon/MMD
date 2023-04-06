function msf_mkdir(folder_path)
% function msf_mkdir(folder_path)
%
% Recursively makes folder_path

if (isempty(folder_path)), return; end

if (isunix)
    
    try
        eval(['!mkdir -p ' folder_path]);
        %Fix to work with spaces in path. DT 20200201
        %Fix stopped working in macOS 10.15.4 Catalina. Reverted to original
        %above without "$". DT 20200626
%         eval(['!mkdir -p $"' folder_path '"']); 
    catch 
        disp('could not make directory');
    end
    
else
    
    [status, ~, ~] = mkdir(folder_path);
    
    if (~status)
        error('Could not make directory');
    end
    
end