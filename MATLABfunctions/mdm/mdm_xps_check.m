function mdm_xps_check(xps)
% function mdm_xps_check(xps)
%
% Check that necessary fields are present in the xps, and that all other
% fields are well formatted

if (~isfield(xps,'n')), error('xps.n is required'); end

f = fieldnames(xps);
for c = 1:numel(f)
    if (strcmp(f{c}, 'n')), continue; end
    if (strcmp(f{c}, 'intent')), continue; end
    if (strcmp(f{c}, 'c_volume')), continue; end
    
    if (size(xps.(f{c}), 1) ~= xps.n)
        error('field %s is of wrong size (%s)', f{c}, num2str(size(xps.(f{c})))); 
    end
    
end