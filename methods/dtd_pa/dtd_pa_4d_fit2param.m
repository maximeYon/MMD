function dps = dtd_pa_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_pa_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
dps = mdm_mfs_load(mfs_fn);

% create parameter maps and save them
sz = size(dps.m);
udtiparam = {'udelta','ufa','ucs','ucl','ucp'};
param = {udtiparam{:},'s0','miso','viso','maniso','vaniso','msaniso','vsaniso'};
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
                    aniso_v = par - perp;
                    saniso_v = aniso_v.^2;

                    miso = iso_v'*w/s0;
                    viso = (iso_v-miso)'.^2*w/s0;
                    maniso = aniso_v'*w/s0;
                    vaniso = (aniso_v-maniso)'.^2*w/s0;
                    msaniso = saniso_v'*w/s0;
                    vsaniso = (saniso_v-msaniso)'.^2*w/s0;

                    dps.s0(ni,nj,nk) = s0;
                    dps.miso(ni,nj,nk) = miso;
                    dps.viso(ni,nj,nk) = viso;
                    dps.maniso(ni,nj,nk) = maniso;
                    dps.vaniso(ni,nj,nk) = vaniso;
                    dps.msaniso(ni,nj,nk) = msaniso;
                    dps.vsaniso(ni,nj,nk) = vsaniso;

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
dps.viso_n = dps.viso./dps.miso.^2;
dps.viso_n(isnan(dps.viso_n)) = 0;
dps.maniso_n = dps.maniso./dps.miso;
dps.maniso_n(isnan(dps.maniso_n)) = 0;
dps.vaniso_n = dps.vaniso./dps.miso.^2;
dps.vaniso_n(isnan(dps.vaniso_n)) = 0;
dps.msaniso_n = dps.msaniso./dps.miso.^2;
dps.msaniso_n(isnan(dps.msaniso_n)) = 0;
dps.vsaniso_n = dps.vsaniso./dps.miso.^4;
dps.vsaniso_n(isnan(dps.vsaniso_n)) = 0;

%2nd moments of P(Deff)
dps.mu2iso = dps.viso;
dps.mu2aniso = 4/5*dps.msaniso;

dps.mdelta = dps.maniso_n/3;

if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end



