function fexi11_check_xps(xps)
% function fexi11_check_xps(xps)
% 
% Checks that all required fields are found in the xps. The required fields
% are:
% 
% mde_b1      - filter b-value
% mde_b2      - detection b-value
% mde_tm12    - mixing time
% s_ind       - series index, yields separate b0's
% mde_b2_ind  - detection b-values, for guessing s0


mdm_xps_check(xps);

required_fields = {...
    'mde_b1', ...   % filter b-value
    'mde_b2', ...   % detection b-value
    'mde_tm12', ... % mixing time from a diffusion point of view
    's_ind', ...    % series index, yields separate s0's
    'mde_b2_ind'};  % unqiue detection b-values, for guessing s0

for c = 1:numel(required_fields)
    if (~isfield(xps, required_fields{c}))
        error('xps.%s not found', required_fields{c});
    end
end
    
