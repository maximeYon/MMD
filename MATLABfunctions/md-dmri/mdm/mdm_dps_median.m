function median_dps = mdm_dps_median(bs_dps)
% function median_dps = mdm_dps_median(bs_dps)
%

ind_bs = find(~cellfun('isempty',bs_dps));
if numel(ind_bs) == 0
    warning('bs_dps is empty')
    return
end

median_dps.nii_h = bs_dps{ind_bs(1)}.nii_h;

sz = ones(1,3);
sz_temp = size(bs_dps{ind_bs(1)}.s0);
sz(1:numel(sz_temp)) = sz_temp;

f = fieldnames(bs_dps{ind_bs(1)});
for c = 1:numel(f)
    if (isstruct(bs_dps{ind_bs(1)}.(f{c})))
        continue;
    elseif iscell(bs_dps{ind_bs(1)}.(f{c}))

        for cbin = 1:numel(bs_dps{ind_bs(1)}.(f{c}))
            fbin = fieldnames(bs_dps{ind_bs(1)}.(f{c}){cbin});
            for cfbin = 1:numel(fbin)
                if (size(bs_dps{ind_bs(1)}.(f{c}){cbin}.(fbin{cfbin}), 1) == sz(1) && ndims(bs_dps{ind_bs(1)}.(f{c}){cbin}.(fbin{cfbin}))<4)                     
                    ptemp = zeros(sz(1),sz(2),sz(3),numel(ind_bs));
                    for nbs = 1:numel(ind_bs)
                        ptemp(:,:,:,nbs) = bs_dps{ind_bs(nbs)}.(f{c}){cbin}.(fbin{cfbin});
                    end
                    median_dps.(f{c}){cbin}.(fbin{cfbin}) = msf_notfinite2zero(nanmedian(ptemp,4));
                end
            end                   
        end

    elseif (size(bs_dps{ind_bs(1)}.(f{c}), 1) == sz(1) && ndims(bs_dps{ind_bs(1)}.(f{c}))<4) 
        ptemp = zeros(sz(1),sz(2),sz(3),numel(ind_bs));
        for nbs = 1:numel(ind_bs)
            ptemp(:,:,:,nbs) = bs_dps{ind_bs(nbs)}.(f{c});
        end
        median_dps.(f{c}) = msf_notfinite2zero(nanmedian(ptemp,4));
    end
end

% Difficult to take median of vectors so the main eigenvector u is
% calculated from the median diffusion tensor elements mdij
if isfield(median_dps,'mdxx')
    % reshape help functions
    sz_reshape  = msf_size(median_dps.s0, 3); 
    g_reshape = @(a,n) reshape(a, prod(sz_reshape(1:3)), n);
    f_reshape = @(a,n) reshape(a, sz_reshape(1), sz_reshape(2), sz_reshape(3), n);
    
    %Voigt format [xx, yy, zz, sqrt(2)*xy, sqrt(2)*xz, sqrt(2)*xz]
    dt_1x6 = cat(4,median_dps.mdxx,median_dps.mdyy,median_dps.mdzz,sqrt(2)*median_dps.mdxy,sqrt(2)*median_dps.mdxz,sqrt(2)*median_dps.mdyz);
    temp_dps = tm_dt_to_dps(g_reshape(dt_1x6, 6)*1e9, median_dps, f_reshape, 0.0001);
    median_dps.u = temp_dps.u;
end

