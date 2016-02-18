function s = mio_mask_simple(s,o)
% function s = mio_mask_simple(s,o)

opt = mdm_opt();

% construct the filename
[~,name] = msf_fileparts(s.nii_fn);
s.mask_fn = fullfile(o, [name '_fit_mask' opt.nii_ext]);

if (exist(s.mask_fn, 'file') && (~opt.do_overwrite))
    return; 
end

[I,h] = mdm_nii_read(s.nii_fn);

% I_max = nanstd(single(abs(I)), [], 4);
% I_max = imfilter(I_max, ones(3,3,3) / 3^3);
% 
% I_std = nanstd(single(I),[],4)

% define the mask from the variation of the signal
I_mean = nanmean(single(abs(I)), 4);
I_mean = imfilter(I_mean, ones(5,5,5) / 5^3);

I_V = zeros(size(I_mean));
for c = 1:size(I,4)
    I_V = I_V + (single(abs(I(:,:,:,c))) - imfilter(single(abs(I(:,:,:,c))), ones(3,3,3) / 3^3)).^2;
end
I_V = I_V / size(I,4);

M = (I_mean ./ sqrt(I_V)) > 2.5; % empirical limit

M = imfilter(single(M), ones(3,3,3)) > 10; % get rid of lonely voxels
M = imfilter(single(M), ones(3,3,3)) > 0;  % fill up the edges

M = mio_mask_fill(M,3);
M = mio_mask_fill(M,2);
M = mio_mask_fill(M,1);
M = mio_mask_fill(M,2);
M = mio_mask_fill(M,3);
M = mio_mask_keep_largest(M);

% write the mask, don't care if we overwrite anything
mdm_nii_write(uint8(M), s.mask_fn, h);