function vasco16_check_xps(xps)
% function vasco16_check_xps(xps)

mdm_xps_check(xps);

if (~isfield(xps, 'alpha2')), error('xps.alpha2 required'); end
if (~isfield(xps, 'b')), error('xps.b required'); end

% This field is required to allow FC and NC data sets to have different
% values of s0
if (~isfield(xps, 's_ind')), error('xps.s_ind required'); end

tmp = unique(xps.s_ind);

if (numel(tmp) ~= 2), error('function built for xps.s_ind) == 1 or 2'); end
if (min(tmp) ~= 1) || (max(tmp) ~= 2)
    error('function built for xps.s_ind) == 1 or 2'); 
end
