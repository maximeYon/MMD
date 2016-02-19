function s = mio_smooth_4d(s, o_path, sigma, opt)
% function s = mio_smooth_4d(s, o_path, sigma, opt)
%
% Smoothes every volume in s.nii_fn and saves this as a new nifti file in
% the folder 'o_path'
%
% The 's' structure is updated with a reference to the new file in s.nii_fn
%
% The smoothing width is controll by 'sigma'

% init 
if (nargin < 4), opt.present = 1; end
opt = mdm_opt(opt);

[~,name] = msf_fileparts(s.nii_fn);
out_nii_fn = fullfile(o_path, [name '_s' opt.nii_ext]);

if (exist(out_nii_fn, 'file') && (~opt.do_overwrite))
    disp('found output, returning');
    s.nii_fn = out_nii_fn;
    return;
end

% Load data
[I,h] = mdm_nii_read(s.nii_fn);

% Create filter
filter_sigma = sigma .* (mean(h.pixdim(2:4)) ./ h.pixdim(2:4));
filter = my_fspecial('gaussian', [0 0 0] + round(max(filter_sigma) * 3)*2+1, ...
    filter_sigma, size(I(:,:,:,1)));

% Determine type of input
tmp = whos('I');
d_type = tmp.class;

% Smooth
for c = 1:size(I,4)
    I_tmp = double(I(:,:,:,c));
    I_tmp = imfilter(I_tmp, filter);
    I(:,:,:,c) = cast(I_tmp, d_type);
end

% save output
s.nii_fn = out_nii_fn;
mdm_nii_write(I, s.nii_fn, h);

end

function filter = my_fspecial(filter_type, filter_size, filter_std, I_size)
% function filter = my_fspecial(filter_type, filter_size, filter_std, I_size)
%
% currently, this function only produce a 3d gaussian filter, but could be
% extended later

if (~strcmp(filter_type, 'gaussian'))
    error('function only defined for gaussian filters');
end

if (numel(I_size) == 2)
    filter = fspecial(filter_type, filter_size(1:2), median(filter_std));
    return;
end

if (numel(filter_size) ~= 3)
    error('only for 3d-filtering for now');
end

x      = zeros(filter_size);
y      = zeros(filter_size);
z      = zeros(filter_size);

for i = 1:filter_size(1)
    for j = 1:filter_size(2)
        for k = 1:filter_size(3)
            x(i,j,k) = i - (filter_size(1)+1)/2;
            y(i,j,k) = j - (filter_size(2)+1)/2;
            z(i,j,k) = k - (filter_size(3)+1)/2;
        end
    end
end

if (numel(filter_std) == 1)
    filter = exp( (- (x.^2 + y.^2 + z.^2) ) / (2 * filter_std^2));
elseif (numel(filter_std) == 3)
    filter = ...
        exp( (- (x.^2 ) ) / (2 * filter_std(1)^2)) .* ...
        exp( (- (y.^2 ) ) / (2 * filter_std(2)^2)) .* ...
        exp( (- (z.^2 ) ) / (2 * filter_std(3)^2));
    
else
    error('strange number of filter stds');
end

filter = filter / sum(filter(:));

end

