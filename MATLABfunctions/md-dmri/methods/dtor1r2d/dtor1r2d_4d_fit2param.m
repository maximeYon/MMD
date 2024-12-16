function dps = dtor1r2d_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtor1r2d_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
opt = dtor1r2d_opt(opt);
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

omega = opt.dtor1r2d.maps_omega;
Nomega = numel(omega);
for nomega = 1:Nomega
    %[dpar,dperp,theta,phi,d0,rpar,rperp,w] = dtor1r2d_4d_m2pars(m);
    [dpar,dperp,theta,phi,r1,r2,w] = dtor1r2d_4d_m2parso(m,omega(nomega));

    sz = size(m);
    nn = size(dpar,4);

    %Calculate derived parameters
    [dxx,dyy,dzz,dxy,dxz,dyz] = dtd_pars2elements(dpar,dperp,theta,phi);
    [diso,daniso,dratio,ddelta,sdaniso,sddelta] = dtd_pars2dpars(dpar,dperp);

    dtor1r2ds = struct('w',w,'dpar',dpar,'dperp',dperp,'theta',theta,'phi',phi,'diso',diso,'daniso',daniso,'ddelta',ddelta,...
        'sdaniso',sdaniso,'sddelta',sddelta,'dratio',dratio,'dxx',dxx,'dyy',dyy,'dzz',dzz,'dxy',dxy,'dxz',dxz,'dyz',dyz,'r1',r1,'r2',r2);

    dps = dtor1r2d_dtor1r2ds2dps(dps, dtor1r2ds);

    % reshape help functions
    sz_reshape  = msf_size(m(:,:,:,1), 3);
    g_reshape = @(a,n) reshape(a, prod(sz_reshape(1:3)), n);
    f_reshape = @(a,n) reshape(a, sz_reshape(1), sz_reshape(2), sz_reshape(3), n);
    dt = cat(4,dps.mdxx,dps.mdyy,dps.mdzz,dps.mdxy,dps.mdxz,dps.mdyz);
    dps = tm_dt_to_dps(g_reshape(dt, 6)*1e9, dps, f_reshape, 0.0001);
%     [min(diso(:)) max(diso(:))]
%     [min(sddelta(:))  max(sddelta(:))]
    %Per-bin statistical measures
    for nbin = 1:numel(opt.dtor1r2d.bin_disomax)    
        if nomega == 1 % Binning done for lowest frequency
            ind_bin = false([sz(1) sz(2) sz(3) nn 10]);
            ind_bin(:,:,:,:,1) = diso >= opt.dtor1r2d.bin_disomin(nbin);
            ind_bin(:,:,:,:,2) = diso <= opt.dtor1r2d.bin_disomax(nbin);
            ind_bin(:,:,:,:,3) = dratio >= opt.dtor1r2d.bin_dratiomin(nbin);
            ind_bin(:,:,:,:,4) = dratio <= opt.dtor1r2d.bin_dratiomax(nbin);
            ind_bin(:,:,:,:,5) = sddelta >= opt.dtor1r2d.bin_sddeltamin(nbin);
            ind_bin(:,:,:,:,6) = sddelta <= opt.dtor1r2d.bin_sddeltamax(nbin);
            ind_bin(:,:,:,:,7) = r1 >= opt.dtor1r2d.bin_r1min(nbin);
            ind_bin(:,:,:,:,8) = r1 <= opt.dtor1r2d.bin_r1max(nbin);
            ind_bin(:,:,:,:,9) = r2 >= opt.dtor1r2d.bin_r2min(nbin);
            ind_bin(:,:,:,:,10) = r2 <= opt.dtor1r2d.bin_r2max(nbin);
    %         ind_bin = false([sz(1) sz(2) sz(3) nn 4]);
    %         ind_bin(:,:,:,:,1) = diso >= opt.dtor1r2d.bin_disomin(nbin);
    %         ind_bin(:,:,:,:,2) = diso <= opt.dtor1r2d.bin_disomax(nbin);
    %         ind_bin(:,:,:,:,3) = sddelta >= opt.dtor1r2d.bin_sddeltamin(nbin);
    %         ind_bin(:,:,:,:,4) = sddelta <= opt.dtor1r2d.bin_sddeltamax(nbin);
            ind_bin_cell{nbin} = all(ind_bin,5);
        end
        ind_bin = ind_bin_cell{nbin};

        dps_bin.no = nbin;
        dtor1r2ds_temp = dtor1r2ds;
        dtor1r2ds_temp.w = dtor1r2ds.w.*ind_bin;
        dps_bin = dtor1r2d_dtor1r2ds2dps(dps_bin, dtor1r2ds_temp);
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

% fields = {'mdiso'; 'msddelta'};
for nbin = 1:numel(opt.dtor1r2d.bin_disomax)    
    for nfield = 1:numel(fields)
        field = fields{nfield};
        parnam = ['d' field 'dnu'];
        dps.bin{nbin}.(parnam) = (dps.omega{end}.bin{nbin}.(field) - dps.omega{1}.bin{nbin}.(field))./(omega(end)-omega(1))*2*pi;
        parnam2 = ['nd' field 'dnu'];
        dps.bin{nbin}.(parnam2) = (dps.omega{end}.bin{nbin}.(field) - dps.omega{1}.bin{nbin}.(field))./(dps.omega{end}.bin{nbin}.(field) + dps.omega{1}.bin{nbin}.(field))*2./(omega(end)-omega(1))*2*pi;
    end
end

% Clean-up to avoid over-filling the memory
% dps
% dps.bin{1}

dps_temp = dps;
clear dps

parnams = {'s0';'mdiso';'msddelta';'mr1';'mr2';'vdiso';'vsddelta';'vr1';'vr2';...
    'cvdisosddelta';'cvdisor1';'cvdisor2';'cvsddeltar1';'cvsddeltar2';'cvr1r2';...
    'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu';...
    'mdxx';'mdyy';'mdzz';'mdxy';'mdxz';'mdyz';'FA'};
for nparnam = 1:numel(parnams)
    parnam = parnams{nparnam};
    dps.(parnam) = single(dps_temp.(parnam));
end

parnams = {'no';'f';'s0';'mdiso';'msddelta';'mr1';'mr2';'vdiso';'vsddelta';'vr1';'vr2';...
    'cvdisosddelta';'cvdisor1';'cvdisor2';'cvsddeltar1';'cvsddeltar2';'cvr1r2';...
    'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu';...
    'mdxx';'mdyy';'mdzz';'mdxy';'mdxz';'mdyz'};
for nbin = 1:numel(opt.dtor1r2d.bin_disomax) 
    for nparnam = 1:numel(parnams)
        parnam = parnams{nparnam};
        dps.bin{nbin}.(parnam) = single(dps_temp.bin{nbin}.(parnam));
    end
end

% dps
% dps.bin{1}
% End clean-up

if (~isempty(dps_fn)) mdm_dps_save(dps, dps.s, dps_fn, opt); end

end

