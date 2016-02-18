function s = mio_powder_average(s, o_path, opt)
% function s = mio_powder_average(s, o_path, opt)
%
% Average over rotations. Image volumes with identical rotations is defined
% from s.xps.a_ind
%
% To do: find a way of keeping track of number of averages per step

% Init
if (nargin < 3), opt.present = 1; end
opt = mdm_opt(opt);

[~,name] = msf_fileparts(s.nii_fn);
out_nii_fn = fullfile(o_path, [name '_pa' opt.nii_ext]);

out_xps_fn = fullfile(o_path, [name '_pa_xps.mat']);

if (exist(out_nii_fn, 'file') && (~opt.do_overwrite))
    disp('found output file, returning');
    load(out_xps_fn);
    s.xps = xps;
    s.nii_fn = out_nii_fn;
    return;
end

% Read data
[I, h] = mdm_nii_read(s.nii_fn);
I = double(I);

% Average image
n = max(s.xps.a_ind);
A = zeros(size(I,1), size(I,2), size(I,3), n);

for c = 1:n
    A(:,:,:,c) = nanmean(I(:,:,:,s.xps.a_ind == c),4);
end

% Average fields in xps
s.xps = rmfield(s.xps, 'n');
f = fieldnames(s.xps);
for i = 1:numel(f)
    for c = 1:n
        xps.(f{i})(c,:) = mean(s.xps.(f{i})(s.xps.a_ind == c, :), 1);
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