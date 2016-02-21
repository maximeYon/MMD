function mfs_fn = mdm_mfs_save(mfs, s, o_fn, opt)

if (nargin < 2), o_fn = []; end

if (isempty(o_fn)), o_fn = fullfile(mdm_tmp_path, 'mfs.mat'); end

% Keep this structure among the model fit parameters in order to be able to
% load nifti headers et c later on
mfs.s = s;

% Save data
msf_mkdir(fileparts(o_fn));
save(o_fn, 'mfs');


mfs_fn = o_fn;
