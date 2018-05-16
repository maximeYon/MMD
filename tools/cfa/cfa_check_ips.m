function cfa_check_ips(ips)
% function cfa_check_ips(ips)


% Checks that all required fields are found in the ips
required_fields = {
    'B0'
    'gm'
    'o'
    'fov'
    'res'
    'ipa'
    'ecs'
    'kpv'
    'r_xyz'
    'T2s'
    };

for c = 1:numel(required_fields)
    if (~isfield(ips, required_fields{c}))
        error('ips.%s not found', required_fields{c});
    end
end
   

% Check that direction vectors are normalized to 1
d_name = {'freq.', 'phase', 'slice'};

for i = 1:3
    n = norm(ips.o(i,:));
    
    if (n<(1-eps)) || (n>(1+eps))
        error('Norm of %s vector (ips.o(%i,:)) is not 1', d_name{i}, i);
    end
end

% Check that r_xyz is congruent with FOV and resolution
[X, Y, Z] = fov2xyz(ips.fov, ips.res); 
r_xyz     = [X(:) Y(:) Z(:)];           % position of all voxel centers [m]

if ~all(size(r_xyz)==size(ips.r_xyz)) || ~all(r_xyz(:) == ips.r_xyz(:))
    warning('Voxel positions are not congruent with FOV and resolution! This may be ok, but check it!')
    % To update r_xyz from FOV and res, try calling ips = cfa_ips_set_xyz_from_fov_res(ips)
end