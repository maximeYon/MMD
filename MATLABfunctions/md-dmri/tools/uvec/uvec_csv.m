clear all;

% Define path for output folders
run_path = fileparts(mfilename('fullpath'));
framework_path = fileparts(fileparts(run_path));
pa_path = fullfile(framework_path,'tools','uvec','repulsion_angles');
out_path = fullfile(framework_path,'tools','uvec','repulsion_csv');

for n = 50

    load(fullfile(pa_path,num2str(n)))

    out_fn = fullfile(out_path,[num2str(n) '.txt']);

    [out_path,out_name,out_ext] = fileparts(out_fn);
    if ~isdir(out_path)
        mkdir(out_path)
    end


    m = [theta'; phi']; 
    dlmwrite(out_fn,m,'precision',3) 
    
end
