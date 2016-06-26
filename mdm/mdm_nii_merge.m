function out_nii_fn = mdm_nii_merge(nii_fn_cell, out_nii_fn)
% function mdm_nii_merge(nii_fn_cell, out_nii_fn)

I = [];
for c = 1:numel(nii_fn_cell)
    
    [I_tmp,h] = mdm_nii_read(nii_fn_cell{c});
    
    if (c == 1)
        I = I_tmp;
    else
        
        if ...
                (size(I,1) ~= size(I_tmp,1)) || ...
                (size(I,2) ~= size(I_tmp,2)) || ...
                (size(I,3) ~= size(I_tmp,3))
            
            error('nii files to be merged are differently sizes');
            
        end
        
        I = cat(4, I, I_tmp);
                
    
    end
end

mdm_nii_write(I, out_nii_fn, h); 