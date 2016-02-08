function mfs_fn = qdti_4d_data2fit(s, o, opt)
% function mfs = qdti_4d_data2fit(s, o, opt)

if (nargin < 3), opt = []; end

opt = qdti_opt(opt);

% Read and reformat data
I  = mdm_nii_read(s.nii_fn);
sz = size(I);
I  = reshape(I, prod(sz(1:3)), sz(4));

% Manage input data; replaces all zero values or below with first value
% above zero in a volume-by-volume approach to save memory
for c = 1:size(I,2)   
    tmp = I(:,c);
    val = min(tmp(tmp(:) > 0));
    tmp(tmp <= 0) = val;
end



% Setup the regression problem: log(I) =   B    *   X
%                               m * n  = (m * 7) (7 * n)
X = [ ones(s.xps.n,1) -s.xps.bt * 1e-9]';

B = log(single(I)) / X;


% Create the model fit structure
mfs.s0 = reshape(real(exp(B(:,1))), sz(1:3));
mfs.dt = reshape(B(:,2:7), [sz(1:3) 6]) * 1e-9;

% Save data
mfs_fn = mdm_mfs_save(mfs, s, o, opt);


    
