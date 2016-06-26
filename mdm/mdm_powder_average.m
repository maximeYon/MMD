function s = mdm_powder_average(s, o_path, opt)
% function s = mdm_powder_average(s, o_path, opt)
%
% Average over rotations. Image volumes with identical rotations is defined
% from s.xps.a_ind
%
% To do: find a way of keeping track of number of averages per step

% Init
if (nargin < 2), o_path = fileparts(s.nii_fn); end
if (nargin < 3), opt.present = 1; end
opt = mdm_opt(opt);
msf_log(['Starting ' mfilename], opt);

if (all(ischar(s))), s = mdm_nii_to_s(s); end

% Build filenames
[~,name] = msf_fileparts(s.nii_fn);
out_nii_fn = fullfile(o_path, [name '_pa' opt.nii_ext]);
out_xps_fn = fullfile(o_path, [name '_pa_xps.mat']);

if (exist(out_nii_fn, 'file') && (~opt.do_overwrite))
    disp(['Skipping, output file already exists: ' out_nii_fn]); 
    s.xps = mdm_xps_load(out_xps_fn); 
    s.nii_fn = out_nii_fn;
    return;
end

% Read data
[I, h] = mdm_nii_read(s.nii_fn);
I = double(I);

% Average image
id = s.xps.a_ind;
if (isfield(s.xps,'s_ind'))
    id = [id s.xps.s_ind];
end
[~,~,id_ind] = unique(id, 'rows');

% get rid of NaNs
tmp = sum(isnan(id),2) > 0;
c_list = unique(id_ind(~tmp)); 
n = numel(c_list);

A = zeros(size(I,1), size(I,2), size(I,3), n);

for c = c_list'
    A(:,:,:,c == c_list) = nanmean(I(:,:,:,id_ind == c),4);
end

if (opt.do_pa_abs), A = abs(A); end

% Average fields in xps
s.xps = rmfield(s.xps, 'n');
f = fieldnames(s.xps);
for i = 1:numel(f)
    for c = c_list'
        
        % allow text fields to be just copied
        if (all(ischar(s.xps.(f{i}))))
            xps.(f{i}) = s.xps.(f{i});
            continue; 
        end
        
        try
            xps.(f{i})(c == c_list,:) = mean(s.xps.(f{i})(id_ind == c, :), 1);
        catch
            if (opt.pa_rethrow_error)
                error('failed powder averaging field %s', f{i});
            else
                warning('failed powder averaging field %s', f{i});
            end
        end
    end
end
s.xps = xps;
s.xps.n = n;

% Build a new file name and save
s.nii_fn = out_nii_fn;
mdm_nii_write(single(A), s.nii_fn, h);

% For good manners, save the xps as well
xps = s.xps;
save(out_xps_fn, 'xps');