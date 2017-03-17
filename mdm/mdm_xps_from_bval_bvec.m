function xps = mdm_xps_from_bval_bvec(bval_fn, bvec_fn, b_delta)
% function xps = mdm_xps_from_bval_bvec(bval_fn, bvec_fn, b_delta)
%
% Create an xps from the b-value and b-vec file resulting from a dcm2nii
% conversion
%
% bval_fn - path to b-value file
% bvec_fn - path to vector file
%
% optional
% b_delta - anisotropy of the b-tensor, range -0.5-1.0 (default 1.0)

if (nargin < 3), b_delta = 1; end

if (nargin == 1) || (isempty(bvec_fn)) % assume it is a nifti filename
    [bval_fn, bvec_fn] = mdm_fn_nii2bvalbvec(bval_fn);
end

if (~exist(bval_fn, 'file')), error('could not find %s', bval_fn); end
if (~exist(bvec_fn, 'file')), error('could not find %s', bvec_fn); end

bval = mdm_txt_read(bval_fn);
xps.b = str2num(bval{1})' * 1e6;
xps.n = numel(xps.b);

if (isempty(bvec_fn)), return; end

bvec = mdm_txt_read(bvec_fn);
assert(numel(bvec) == 3, 'strange bvec file');

gdir = [str2num(bvec{1}); str2num(bvec{2}); str2num(bvec{3})];

% compute b-tensors from b-values, b_delta value(s) and symmetry axis
bt  = tm_tpars_to_1x6(xps.b, b_delta, gdir');
xps = mdm_xps_from_bt(bt);

% store the direction as well
xps.u = gdir';


