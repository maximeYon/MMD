function my_function_python_denoising_2D(data_path)

if ispc == 1
    %% suppose Anaconda3 installation
    pyExec = 'C:\Users\User\anaconda3\';
    pyRoot = fileparts(pyExec);
    p = getenv('PATH');
    p = strsplit(p, ';');
    addToPath = {
        pyRoot
        fullfile(pyRoot, 'Library', 'mingw-w64', 'bin')
        fullfile(pyRoot, 'Library', 'usr', 'bin')
        fullfile(pyRoot, 'Library', 'bin')
        fullfile(pyRoot, 'Scripts')
        fullfile(pyRoot, 'bin')
        };
    p = [addToPath(:); p(:)];
    p = unique(p, 'stable');
    p = strjoin(p, ';');
    setenv('PATH', p);
    
    system('conda activate base');
    system(['python C:\Users\User\Mon_Drive\Python\DESIGNER_denoising\my_function_apply_den_2D.py ' data_path]);
else
    system(['python C:\python_path\my_function_apply_den_2D.py ' data_path]);
end
end