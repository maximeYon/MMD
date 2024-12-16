function xps = mdm_xps_calc_btpars(xps)
% function xps = mdm_xps_calc_btpars(xps)

for c = 1:xps.n
    t = tm_1x6_to_3x3(xps.bt(c,:));
    
    t = tm_3x3_to_tpars(t);
    
    f = fieldnames(t);
    for c_field = 1:numel(f)
        xps.(['b_' f{c_field}])(c,:) = t.(f{c_field});
    end
    
end