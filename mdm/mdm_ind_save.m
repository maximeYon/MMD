function ind_fn = mdm_ind_save(ind, ind_fn)
% function ind_fn = mdm_ind_save(ind, ind_fn)
%
% Saves the model fit structure

% Save data
msf_mkdir(fileparts(ind_fn));
save(ind_fn, 'ind');

