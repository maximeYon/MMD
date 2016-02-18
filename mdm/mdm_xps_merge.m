function xps = mdm_xps_merge(xps_cell)
% function xps = mdm_xps_merge(xps_cell)


for i = 1:numel(xps_cell)
    if (isfield(xps_cell{i}, 's_ind'))
        error('merging xps structures with existing s_ind can be ambigious');
    end
    xps_cell{i}.s_ind = zeros(xps_cell{i}.n, 1) + i;
end

xps = xps_cell{1};

f = fieldnames(xps_cell{1});

for i = 2:numel(xps_cell)
    
    f2 = fieldnames(xps_cell{i});
    
    if (numel(f) ~= numel(f2))
        error('different number of fields present in the xps to be merged');
    end
    
    for c = 1:numel(f)
        
        if (~strcmp(f{c}, f2{c}))
            error('different fields present in the xps to be merged');
        end
            
        xps.(f{c}) = cat(1, xps.(f{c}), xps_cell{i}.(f{c}));
        
        
    end
    
end


xps.n = sum(xps.n);