function mdm_mkdir(folder_path)

if (isunix)
    
    try
        eval(['!mkdir -p ' folder_path]);
    catch 
        disp('could not make directory');
    end
    
else
    
    [status, ~, ~] = mkdir(folder_path);
    
    if (~status)
        error('Could not make directory');
    end
    
end