function fn = dti_nls_4d_fit2param(mfs_fn, o_path, opt)
% function fn = dti_nls_4d_fit2param(mfs_fn, o_path, opt)
%
% Creates meaningful parameters from the model fit structure
%
% Input:
%
% mfs_fn - Path to .mat file with model fit structure
%
% o_path - Output path, where parameters maps are stored
%
% opt    - Options, optional argument
%
%
% Output:
%
% fn     - A cell array with paths to parameter maps that were written

if (nargin < 3), opt = []; end
    
opt = mdm_opt(opt);

mfs = mdm_mfs_load(mfs_fn);
    
mfs.s0(:,:,:,1)   = mfs.m(:,:,:,1); 
mfs.dt(:,:,:,2:7) = mfs.m(:,:,:,2:7); 

% for now, use qdti's algorithm
fn = qdti_4d_fit2param(mfs, o_path, opt);


end