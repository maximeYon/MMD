function mfs_fn = mdm_mfs_save(mfs, s, mfs_fn, opt)
% function mfs_fn = mdm_mfs_save(mfs, s, mfs_fn, opt)
%
% Saves the model fit structure

if (nargin < 2), mfs_fn = []; end
if (nargin < 3), opt.present = 1; end

if (isempty(mfs_fn)), mfs_fn = fullfile(mdm_tmp_path, 'mfs.mat'); end

% Keep this structure among the model fit parameters in order to be able to
% load nifti headers et c later on
mfs.s = s;

% Save data
msf_mkdir(fileparts(mfs_fn));

s = whos('mfs');
sizeGB = s.bytes / (1024^3);
if sizeGB >= 2
    save(mfs_fn, 'mfs','-v7.3');
else
    save(mfs_fn, 'mfs');
end


