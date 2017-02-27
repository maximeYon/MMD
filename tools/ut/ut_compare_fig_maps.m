function ut_compare_fig_maps(a, b, k, method, do_show)
% function ut_compare_fig_maps(a, b, k, method)
%
% Compare image volumes in folders 'a' and 'b', in transversal slice 'k'
% 
% The maps that are obtained are fetched from opt.(method).fig_maps.
%
% Optional
% do_show - show actual image maps (with pause, disables errors)
%           


if (nargin < 5), do_show = 0; end

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
    
    A(isnan(A(:))) = 0;
    B(isnan(B(:))) = 0;
    
    % compare the two: allow a difference of 1% of max image value in at
    % most 1% of all voxels
    n = sum( abs(A(:) - B(:)) > 1e-2 * max(q(A),q(B)));
    n_tol = 0.01 * sum( (A(:) ~= 0) & (B(:) ~= 0) );
    if (n > n_tol)
        
        if (do_show)
            f = @(varargin) warning(varargin{:});
        else
            f = @(varargin) error(varargin{:});
        end
        
        f('%s, %s, analysed file differs (%s)', mfilename, method, map);
    end        
    
    
    if (do_show)
        
        c_min = min(min(A(:)), min(B(:)));
        c_min = min(0, c_min);
        c_max = max(quantile(A(:), 0.99), quantile(B(:), 0.99));
        
        subplot(2,2,1);
        msf_imagesc(A);
        caxis([c_min c_max]);
        colorbar;
        title(strrep(p_fn, '_', '\_'));
        
        subplot(2,2,2);
        msf_imagesc(B);
        caxis([c_min c_max]);
        colorbar;
        
        subplot(2,2,3);
        msf_imagesc(A-B);
        caxis([-1 1] * max(abs([c_min c_max])) * 0.1);
        colorbar;
        
        pause;
    end
        
        
    
end
