function M = mdm_mask_load(s, opt)
% function M = mdm_mask_load(s, opt)
%
% Load the mask from s.mask_fn
%
% Optionally, mask using i_range, j_range and k_range in opt

h = mdm_nii_read_header(s.nii_fn);

sz = h.dim(2:4)';

if (~isfield(s, 'mask_fn'))
    M = ones(sz); 
else
    M = mdm_nii_read(s.mask_fn);
    
    if (~all(size(M) == sz))
        error(['Size of mask different from data (' num2str(size(M)) ' vs ' num2str(sz)]);
    end
end

% Optional masking of i,j,k ranges
M2 = ones(sz);
if (isfield(opt, 'i_range')), M3 = zeros(sz); M3(opt.i_range,:,:) = 1; M2 = M2 & M3; end
if (isfield(opt, 'j_range')), M3 = zeros(sz); M3(:,opt.j_range,:) = 1; M2 = M2 & M3; end
if (isfield(opt, 'k_range')), M3 = zeros(sz); M3(:,:,opt.k_range) = 1; M2 = M2 & M3; end

M = M & M2;

    
    

