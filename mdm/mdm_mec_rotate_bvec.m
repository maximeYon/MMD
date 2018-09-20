function xps = mdm_mec_rotate_bvec(xps, tpm_fn, p_fn)
% function mdm_mec_rotate_bvec(xps, tpm_fn, p_fn)

% First check that affine dti transform was used
p = elastix_p_read(p_fn);

if (~strcmp(p.Transform, '"AffineDTITransform"'))
    error('Expected AffineDTITransform to be used');
end

% Get the transform parameters
tpm = elastix_tpm_read(tpm_fn);


% Prepare for matrix rotation
for c = 1:xps.n
    tmp = reshape(tpm(:,c), 3, 4)';
    tmp(2,:) = 0; % translation
    tmp(3,:) = 1; % scales
    tmp(4,:) = 0; % shears
    tmp_tmat = elastix_param2tmat(tmp);
    tmp_tmat = inv(tmp_tmat); % go back to original space
    tmp_M = tmp_tmat(1:3,1:3); % removes translations
    
    % tmp * tmp' = eye, unless there are problem
    R = ((tmp_M * tmp_M')^(-1/2)) * tmp_M;
    
    % Apply
    if (isfield(xps, 'u'))
        xps.u(c,:) = (R * xps.u(c,:)')';
    end
    
    if (isfield(xps, 'bt'))
        xps.bt(c,:) = tm_3x3_to_1x6( R * tm_1x6_to_3x3(xps.bt(c,:)) * R' );
    end
    
    
    
end

