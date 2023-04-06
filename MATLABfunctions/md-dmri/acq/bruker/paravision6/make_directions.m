% Converts repulsion directions to text string for pasting into PV methods
% file in the fields DwDir and DwDirInput

clear all;

% Define path for output folders
run_path = pwd;
framework_path = fileparts(fileparts(fileparts(run_path)));
pa_path = fullfile(framework_path,'tools','uvec','repulsion_angles');
out_path = fullfile(run_path, 'directions');
msf_mkdir(out_path);

for n = 10:128

    load(fullfile(pa_path,num2str(n)))

    out_fn = fullfile(out_path,[num2str(n) '.txt']);

    [out_path,out_name,out_ext] = fileparts(out_fn);
    if ~isdir(out_path)
        mkdir(out_path)
    end

    fid = fopen(out_fn,'w');

    formatspec = '%8.15f ';
    for nangle = 1:numel(theta)
        out_mat = [cos(phi(nangle)).*sin(theta(nangle)) sin(phi(nangle)).*sin(theta(nangle)) cos(theta(nangle))];
        out_mat = out_mat./repmat(sqrt(out_mat(1)^2+out_mat(2)^2+out_mat(3)^2),[1 3]);
        for nxyz = 1:3
            fprintf(fid, formatspec, out_mat(nxyz));
        end
    end
    fclose(fid);

end
