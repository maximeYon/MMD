function ut_compare_fig_maps(a, b, k, method)
% function ut_compare_fig_maps(a, b, k, method)

opt = feval([method '_opt']); 

q = @(x) quantile(x( x(:) ~= 0), 0.95);

for c = 1:numel(opt.(method).fig_maps)
    
    map = opt.(method).fig_maps{c};
    p_fn = [opt.(method).fig_prefix '_' map '.nii.gz'];
    
    if (~exist(fullfile(a, p_fn), 'file'))
        error('%s, %s output file not found (%s)', mfilename, method, map);
    end

    if (~exist(fullfile(b, p_fn), 'file'))
        error('%s, %s true output file not found (%s)', mfilename, method, map);
    end
    
    A = mdm_nii_read(fullfile(a, p_fn));
    B = mdm_nii_read(fullfile(b, p_fn));
    
    % select the slice of interest
    A = A(:,:,k);
    B = B(:,:,k);
    
    % compare the two
    if (any( abs(A(:) - B(:)) > 1e-3 * max(q(A),q(B))))
        error('%s, gamma analysed file differs (%s)', mfilename, map);
    end
        
    
end
