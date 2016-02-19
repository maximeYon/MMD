function fexi11_check_xps(xps)
% function fexi11_check_xps(xps)
% 
% checks that all required fields are found in the xps


required_fields = {...
    'mde_b1', 'mde_b2', 'mde_tm12', 's_ind', 'mde_b2_ind'};

for c = 1:numel(required_fields)
    if (~isfield(xps, required_fields{c}))
        error('xps.%s not found', required_fields{c});
    end
end
    
