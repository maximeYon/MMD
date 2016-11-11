function xps = mdm_xps_from_bval_bvec(bval_fn, bvec_fn, b_delta)
% function xps = mdm_xps_from_bval_bvec(bval_fn, bvec_fn)

if (nargin < 3), b_delta = 1; end

if (nargin == 1) || (isempty(bvec_fn)) % assume it is a nifti filename
    [nii_path, nii_name] = msf_fileparts(bval_fn);
    bval_fn = fullfile(nii_path, [nii_name '.bval']);
    bvec_fn = fullfile(nii_path, [nii_name '.bvec']);
end

if (~exist(bval_fn, 'file')), error('could not find %s', bval_fn); end

bval = mdm_txt_read(bval_fn);
xps.b = str2num(bval{1})' * 1e6;
xps.n = numel(xps.b);


if (isempty(bvec_fn))
    return;
end

bvec = mdm_txt_read(bvec_fn);
assert(numel(bvec) == 3, 'strange bvec file');

gdir = [str2num(bvec{1}); str2num(bvec{2}); str2num(bvec{3})];

xps.bt = tm_1x3_to_1x6(xps.b, zeros(size(xps.b)), gdir');


if (b_delta == 0)
    xps.bt = repmat([1/3 1/3 1/3 0 0 0], size(xps.bt,1), 1) .* repmat(xps.b, 1, 6);
elseif (b_delta ~= 1)
    error('not yet implemented');
end

xps = mdm_xps_from_bt(xps.bt);

xps.u = gdir';


