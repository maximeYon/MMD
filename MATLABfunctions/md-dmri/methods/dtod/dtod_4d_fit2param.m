function dps = dtod_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtod_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
opt = dtod_opt(opt);
%dps = mdm_mfs_load(mfs_fn);

% Hack to allow mgui to access this function
if ischar(mfs_fn)
    dps = mdm_mfs_load(mfs_fn);
else
    m = mfs_fn;
    dps.m = m;
end

% create parameter maps and save them

m = dps.m;
dps = rmfield(dps,'m');

omega = opt.dtod.maps_omega;
Nomega = numel(omega);
for nomega = 1:Nomega
    %[dpar,dperp,theta,phi,d0,rpar,rperp,w] = dtod_4d_m2pars(m);
    [dpar,dperp,theta,phi,w] = dtod_4d_m2parso(m,omega(nomega));

    sz = size(m);
    nn = size(dpar,4);

    %Calculate derived parameters
    [dxx,dyy,dzz,dxy,dxz,dyz] = dtd_pars2elements(dpar,dperp,theta,phi);
    [diso,daniso,dratio,ddelta,sdaniso,sddelta] = dtd_pars2dpars(dpar,dperp);

    dtods = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,'diso',diso,'daniso',daniso,'ddelta',ddelta,...
        'sdaniso',sdaniso,'sddelta',sddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,'dxy',dxy,'dxz',dxz,'dyz',dyz);

    dps = dtod_dtods2dps(dps, dtods);

    % reshape help functions
    sz_reshape  = msf_size(m(:,:,:,1), 3);
    g_reshape = @(a,n) reshape(a, prod(sz_reshape(1:3)), n);
    f_reshape = @(a,n) reshape(a, sz_reshape(1), sz_reshape(2), sz_reshape(3), n);
    dt = cat(4,dps.mdxx,dps.mdyy,dps.mdzz,dps.mdxy,dps.mdxz,dps.mdyz);
    dps = tm_dt_to_dps(g_reshape(dt, 6)*1e9, dps, f_reshape, 0.0001);

    %Per-bin statistical measures
    for nbin = 1:numel(opt.dtod.bin_disomax) 
        if nomega == 1
            ind_bin = false([sz(1) sz(2) sz(3) nn 6]);
        %     ind_bin = false([sz(1) sz(2) sz(3) nn 8]);
            ind_bin(:,:,:,:,1) = diso >= opt.dtod.bin_disomin(nbin);
            ind_bin(:,:,:,:,2) = diso <= opt.dtod.bin_disomax(nbin);
            ind_bin(:,:,:,:,3) = dratio >= opt.dtod.bin_dratiomin(nbin);
            ind_bin(:,:,:,:,4) = dratio <= opt.dtod.bin_dratiomax(nbin);
            ind_bin(:,:,:,:,5) = sddelta >= opt.dtod.bin_sddeltamin(nbin);
            ind_bin(:,:,:,:,6) = sddelta <= opt.dtod.bin_sddeltamax(nbin);
        %     ind_bin(:,:,:,:,7) = r1 >= opt.dtod.bin_r1min(nbin);
        %     ind_bin(:,:,:,:,8) = r1 <= opt.dtod.bin_r1max(nbin);
            ind_bin_cell{nbin} = all(ind_bin,5);
        end
        ind_bin = ind_bin_cell{nbin};

        dps_bin.no = nbin;
        dtods_temp = dtods;
        dtods_temp.w = dtods.w.*ind_bin;
        dps_bin = dtod_dtods2dps(dps_bin, dtods_temp);
        dps_bin.f = dps_bin.s0./dps.s0;
        dps.bin{nbin} = dps_bin;
    end

    dps.omega{nomega} = dps;
    if nomega > 1
        dps.omega{nomega} = rmfield(dps.omega{nomega},'omega');
    end

end

dps_temp = dps;
dps = dps_temp.omega{1}; % Lowest frequency
dps.omega = dps_temp.omega;

%Rate of change with nu = omega/(2*pi);
fields = {'mdiso'; 'msddelta'; 'vdiso'; 'vsddelta'; 'cvdisosddelta'};
for nfield = 1:numel(fields)
    field = fields{nfield};
    parnam = ['d' field 'dnu'];
    dps.(parnam) = (dps.omega{end}.(field) - dps.omega{1}.(field))./(omega(end)-omega(1))*2*pi;
    parnam2 = ['nd' field 'dnu'];
    dps.(parnam2) = (dps.omega{end}.(field) - dps.omega{1}.(field))./(dps.omega{end}.(field) + dps.omega{1}.(field))*2./(omega(end)-omega(1))*2*pi;
end

%fields = {'mdiso'; 'msddelta'};
for nbin = 1:numel(opt.dtod.bin_disomax)    
    for nfield = 1:numel(fields)
        field = fields{nfield};
        parnam = ['d' field 'dnu'];
        dps.bin{nbin}.(parnam) = (dps.omega{end}.bin{nbin}.(field) - dps.omega{1}.bin{nbin}.(field))./(omega(end)-omega(1))*2*pi;
        parnam2 = ['nd' field 'dnu'];
        dps.bin{nbin}.(parnam2) = (dps.omega{end}.bin{nbin}.(field) - dps.omega{1}.bin{nbin}.(field))./(dps.omega{end}.bin{nbin}.(field) + dps.omega{1}.bin{nbin}.(field))*2./(omega(end)-omega(1))*2*pi;
    end
end

if (~isempty(dps_fn)) mdm_dps_save(dps, dps.s, dps_fn, opt); end

end

