function s = mio_mask_mic(s, o, opt)
% function s = mio_mask_mic(s, o, maskopt)

if (nargin < 3)
    opt.mask.b0_ind = 1;
    opt.mask.thresh = .1;
end
%opt = mdm_opt();

[I,h] = mdm_nii_read(s.nii_fn);

% % define the mask from the variation of the signal
% cv = nanstd(single(I),[],4) ./ nanmean(single(I), 4);
% cvlim = 0.8;
% M = cv > cvlim; % empirical limit
% %figure(1), clf, imagesc(cv'), set(gca,'YDir','normal','Clim',[0 cvlim]), axis equal, axis tight, axis off, pause, return

% define the mask from 10% threshold of b0 image
I0 = I(:,:,:,opt.mask.b0_ind);
Imax = max(reshape(I0,[numel(I0) 1]));
M = I0/Imax > opt.mask.thresh;
%figure(1), clf, imagesc(I0'/Imax), set(gca,'YDir','normal','Clim',[0 opt.mask.thresh]), axis equal, axis tight, axis off, pause, return


% M = imfilter(single(M), ones(3,3,3)) > 10; % get rid of lonely voxels
% M = imfilter(single(M), ones(3,3,3)) > 0;  % fill up the edges

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