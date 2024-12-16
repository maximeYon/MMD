function dtd_check_xps(xps, opt)
% function dtd_check_xps(xps, opt)
% 
% checks that all required fields are found in the xps

mdm_xps_check(xps);

required_fields = {...
    'bt'};

for c = 1:numel(required_fields)
    if (~isfield(xps, required_fields{c}))
        error('xps.%s not found', required_fields{c});
    end
end
    
