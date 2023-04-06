function dps = dtd_pa_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_pa_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
%dps = mdm_mfs_load(mfs_fn);

% Hack to allow mgui to access this function
if ischar(mfs_fn)
    dps = mdm_mfs_load(mfs_fn);
else
    m = mfs_fn;
    dps.m = m; dps.nii_h = []; dps.mask = ones(size(m,1),size(m,2),size(m,3));
end

% create parameter maps and save them
sz = size(dps.m);
udtiparam = {'udelta','ucs','ucp','ucl','ufa'};
param = {udtiparam{:},'s0','mdiso','vdiso','mdaniso','vdaniso','msdaniso','vsdaniso'};
for nparam = 1:numel(param)
    eval(['dps.' param{nparam} ' = zeros([sz(1) sz(2) sz(3)]);']);
end

for nk = 1:sz(3)
    for nj = 1:sz(2)
        for ni = 1:sz(1)
            if dps.mask(ni,nj,nk)
                m = squeeze(dps.m(ni,nj,nk,:))';
                dtd = dtd_pa_m2dtd(m);
                [n,par,perp,w] = dtd_pa_dist2par(dtd);
                                
                if n > 0
                    s0 = ones(1,n)*w;
                    iso_v = (par + 2*perp)/3;
                    aniso_v = (par - perp)/3;
                    saniso_v = aniso_v.^2;

                    mdiso = iso_v'*w/s0;
                    vdiso = (iso_v-mdiso)'.^2*w/s0;
                    mdaniso = aniso_v'*w/s0;
                    vdaniso = (aniso_v-mdaniso)'.^2*w/s0;
                    msdaniso = saniso_v'*w/s0;
                    vsdaniso = (saniso_v-msdaniso)'.^2*w/s0;

                    dps.s0(ni,nj,nk) = s0;
                    dps.mdiso(ni,nj,nk) = mdiso;
                    dps.vdiso(ni,nj,nk) = vdiso;
                    dps.mdaniso(ni,nj,nk) = mdaniso;
                    dps.vdaniso(ni,nj,nk) = vdaniso;
                    dps.msdaniso(ni,nj,nk) = msdaniso;
                    dps.vsdaniso(ni,nj,nk) = vsdaniso;

                    udtd = [par'; perp'; zeros(n,1)'; zeros(n,1)'; w'];
                    udtd = [n; udtd(:)];
                    [udtd_nx6,w] = dtd_dist2nx6w(udtd);
                    udt1x6 = (udtd_nx6'*w)'/s0;
                    udt3x3 = tm_1x6_to_3x3(udt1x6);
                    udt = tm_3x3_to_tpars(udt3x3);

                    dps.udelta(ni,nj,nk) = udt.delta;
                    dps.ufa(ni,nj,nk) = udt.fa;
                    dps.ucs(ni,nj,nk) = udt.cs;
                    dps.ucl(ni,nj,nk) = udt.cl;
                    dps.ucp(ni,nj,nk) = udt.cp;
                end                
            end
        end
    end
end

%normalized parameters
dps.nvdiso = dps.vdiso./dps.mdiso.^2;
dps.nvdiso(isnan(dps.nvdiso)) = 0;
dps.nmdaniso = dps.mdaniso./dps.mdiso;
dps.nmdaniso(isnan(dps.nmdaniso)) = 0;
dps.nvdaniso = dps.vdaniso./dps.mdiso.^2;
dps.nvdaniso(isnan(dps.nvdaniso)) = 0;
dps.nmsdaniso = dps.msdaniso./dps.mdiso.^2;
dps.nmsdaniso(isnan(dps.nmsdaniso)) = 0;
dps.nvsdaniso = dps.vsdaniso./dps.mdiso.^4;
dps.nvsdaniso(isnan(dps.nvsdaniso)) = 0;

%2nd moments of P(Deff)
dps.mu2iso = dps.vdiso;
dps.mu2aniso = 4/5*dps.msdaniso;

dps.mdelta = dps.nmdaniso/3;

dps.mdiso = dps.mdiso;
dps.vdiso = dps.vdiso;
dps.msdaniso = dps.msdaniso;
dps.nvdiso = dps.nvdiso;
dps.nmsdaniso = dps.nmsdaniso;

dps.Vl = 5/2 * 4/5*dps.msdaniso;
dps.MKi = 3 * dps.nvdiso; % Multiply by 3 to get kurtosis
dps.MKa = 3 * 4/5*dps.nmsdaniso;

% Calculate uFA. Take real component to avoid complex values due to
% sqrt of negative variances.
dps.ufa_old = real(sqrt(3/2) * sqrt(1./(dps.mdiso.^2./dps.Vl+1))); % Lasic (2014)
dps.ufa     = real(sqrt(3/2) * sqrt( dps.Vl ./ (dps.Vl + dps.vdiso + dps.mdiso.^2) )); % Szczepankiewicz (2016)

dps.signaniso = sign(dps.nmdaniso);

if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end



