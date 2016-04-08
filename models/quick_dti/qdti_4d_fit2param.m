function fn = qdti_4d_fit2param(mfs_fn, o_path, opt)
% function fn = qdti_4d_fit2param(mfs_fn, o_path, opt)
%
% 2do: Update this function to agree with structure of dti_euler

if (nargin < 3), opt.present = 1; end

% Load data & setup option structure
mfs = mdm_mfs_load(mfs_fn);
h   = mfs.nii_h; 
opt = mdm_opt(opt);

m2v = @(x)   reshape(x, size(x,1) * size(x,2) * size(x,3), size(x,4));
v2m = @(x,r) reshape(x, r(1), r(2), r(3), numel(x) / prod(r(1:3)));

% create parameter maps and save them
n_map = 3;
fn = cell(1,n_map);
for c = 1:n_map
    
    min_max = [-inf inf];
    switch (c)
        
        case 1 % fractional anisotropy
            param = 'fa';
            
            [~,E_shear] = dtd_6x6_iso();
            
            md = dtd_inner(m2v(mfs.dt), dtd_3x3_to_1x6(dtd_3x3_iso()));
            vl = dtd_inner(dtd_1x6_to_1x21(m2v(mfs.dt)), dtd_6x6_to_1x21(E_shear));
            
            x = v2m(sqrt( (3/2) * (1 + md.^2 ./ vl ).^(-1)), size(mfs.dt));
            clear md vl;
            
            min_max = [0 1];
            
        case 2 % mean diffusivity
            param = 'md';
            
            x = v2m(dtd_inner(m2v(mfs.dt), dtd_3x3_to_1x6(dtd_3x3_iso())), size(mfs.dt));
            x = x * 1e9;
            
            min_max = [0 4];
            
        case 3 % s0 
            param = 's0';
            
            x = mfs.s0;
            min_max = [0 inf];
    end
    
    fn{c} = fullfile(o_path, ['qdti_' param opt.nii_ext]);
    
    % make sure the min_max field is there
    opt.qdti.(param).present = 1;
    opt.qdti.(param) = msf_ensure_field(opt.qdti.(param), 'min_max', min_max);
    
    % cut values above/below min/max limits
    x = mio_min_max_cut(x, opt.qdti.(param).min_max);
    
    % write file
    mdm_nii_write(x, fn{c}, h);
end




