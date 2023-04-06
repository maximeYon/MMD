function dti_nls_check_xps(xps)
% function dti_nls_check_xps(xps)
% 
% checks that all required fields are found in the xps

if (~isfield(xps, 'bt')), error('xps.bt not found'); end

