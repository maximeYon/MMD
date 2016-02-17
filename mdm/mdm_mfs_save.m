function mfs_fn = mdm_mfs_save(mfs, s, o, opt)

if (nargin < 2), o = []; end

if (isempty(o)), o = fullfile(mdm_tmp_path, 'mfs.mat'); end

% Keep this structure among the model fit parameters in order to be able to
% load nifti headers et c later on
mfs.s = s;

% Save data
msf_mkdir(fileparts(o));
save(o, 'mfs');


mfs_fn = o;
