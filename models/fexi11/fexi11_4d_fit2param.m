function fn = fexi_nls_4d_fit2param(mfs_fn, o_path, opt)
% function fn = fexi_nls_4d_fit2param(mfs_fn, o_path, opt)
%
% Creates meaningful parameters from the model fit structure
%
% Input:
%
% mfs_fn - Path to .mat file with model fit structure
%
% o_path - Output path, where parameters maps are stored
%
% opt    - Options, optional argument
%
%
% Output:
%
% fn     - A cell array with paths to parameter maps that were written

if (nargin < 3), opt = []; end
    

  
opt = mdm_opt(opt);
mfs = mdm_mfs_load(mfs_fn);
h   = mfs.nii_h; 

% create parameter maps and save them
n_map = 3;
fn = cell(1,n_map);
for c = 1:n_map
    
    min_max = [-inf inf];
    switch (c)
        
        case 1
            param = 'ADC';
            min_max = [0 4];
            x = mfs.m(:,:,:,c)* 1e9;

        case 2
            param = 'sigma';
            min_max = [0 1];
            x = mfs.m(:,:,:,c);
            
        case 3
            param = 'AXR';
            min_max = [0 20];
            x = mfs.m(:,:,:,c);
            

    end
    
    fn{c} = fullfile(o_path, ['fexi16_p_' param opt.nii_ext]);
    
    % make sure the min_max field is there
    opt.fexi16.(param).present = 1;
    opt.fexi16.(param) = msf_ensure_field(opt.fexi16.(param), 'min_max', min_max);
    
    % cut values above/below min/max limits
    x = mio_min_max_cut(x, opt.fexi16.(param).min_max);
    
    % write file
    mdm_nii_write(x, fn{c}, h);
end
