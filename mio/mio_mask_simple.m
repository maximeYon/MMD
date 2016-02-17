function s = mio_mask_simple(s,o)
% function s = mio_mask_simple(s,o)

opt = mdm_opt();

[I,h] = mdm_nii_read(s.nii_fn);

% define the mask from the variation of the signal
cv = nanstd(single(I),[],4) ./ nanmean(single(I), 4);

M = cv > 0.25; % empirical limit

M = imfilter(single(M), ones(3,3,3)) > 10; % get rid of lonely voxels
M = imfilter(single(M), ones(3,3,3)) > 0;  % fill up the edges

M = mio_mask_fill(M,3);
M = mio_mask_fill(M,2);
M = mio_mask_fill(M,1);
M = mio_mask_fill(M,2);
M = mio_mask_fill(M,3);
M = mio_mask_keep_largest(M);

% construct the filename
[~,name] = fileparts(s.nii_fn);

s.mask_fn = fullfile(o, [name '_fit_mask' opt.nii_ext]);

% write the mask, don't care if we overwrite anything
mdm_nii_write(uint8(M), s.mask_fn, h);