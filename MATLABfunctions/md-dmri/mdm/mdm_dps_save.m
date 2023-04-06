function dps_fn = mdm_dps_save(dps, s, dps_fn, opt)
% function dps_fn = mdm_dps_save(dps, s, dps_fn, opt)
% 
% Saves derived parameter structure

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt.present = 1; end

if (isempty(dps_fn)), dps_fn = fullfile(mdm_tmp_path, 'dps.mat'); end

% Keep this structure among the model fit parameters in order to be able to
% load nifti headers et c later on
dps.s = s;

% Save data
msf_mkdir(fileparts(dps_fn));
save(dps_fn, 'dps');


