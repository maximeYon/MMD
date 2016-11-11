function dps = dtd_pa_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_pa_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
dps = mdm_mfs_load(mfs_fn);

% create parameter maps and save them
sz = size(dps.m);
param = {'s0','miso','viso','maniso','vaniso','mvlambda','vvlambda'};
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
                    miso = iso_v'*w/s0;
                    viso = (iso_v-miso)'.^2*w/s0;
                    aniso_v = par - perp;
                    maniso = aniso_v'*w/s0;
                    vaniso = (aniso_v-maniso)'.^2*w/s0;
                    vlambda_v = 2/9*aniso_v.^2;
                    mvlambda = vlambda_v'*w/s0;
                    vvlambda = (vlambda_v-mvlambda)'.^2*w/s0;

                    dps.s0(ni,nj,nk) = s0;
                    dps.miso(ni,nj,nk) = miso;
                    dps.viso(ni,nj,nk) = viso;
                    dps.maniso(ni,nj,nk) = maniso;
                    dps.vaniso(ni,nj,nk) = vaniso;
                    dps.mvlambda(ni,nj,nk) = mvlambda;
                    dps.vvlambda(ni,nj,nk) = vvlambda;
                end                
            end
        end
    end
end

dps.mu2iso = dps.viso;
dps.mu2aniso = 2/5*dps.mvlambda;
dps.mu2aniso(isnan(dps.mu2aniso)) = 0;
dps.mu2tot = dps.mu2iso + dps.mu2aniso;

dps.kiso = dps.mu2iso./dps.miso.^2;
dps.kiso(isnan(dps.kiso)) = 0;
dps.kaniso = dps.mu2aniso./dps.miso.^2;
dps.kaniso(isnan(dps.kaniso)) = 0;
dps.ktot = dps.mu2tot./dps.miso.^2;
dps.ktot(isnan(dps.ktot)) = 0;

dps.ciso = dps.viso./dps.miso.^2;
dps.ciso(isnan(dps.ciso)) = 0;
dps.mdelta = dps.maniso./(3*dps.miso);
dps.mdelta(isnan(dps.mdelta)) = 0;
dps.vdelta = dps.vaniso./(3*dps.miso).^2;
dps.vdelta(isnan(dps.vdelta)) = 0;

dps.cmvlambda = dps.mvlambda./dps.miso.^2;
dps.cmvlambda(isnan(dps.cmvlambda)) = 0;
dps.ufa = sqrt(3/2)*sqrt(1./(1./dps.cmvlambda+1));
dps.ufa(isnan(dps.ufa)) = 0;
dps.cmu = dps.ufa.^2;

dps.cvvlambda = dps.vvlambda./dps.miso.^4;
dps.cvvlambda(isnan(dps.cvvlambda)) = 0;

if (~isempty(dps_fn)), mdm_dps_save(dps, dps.s, dps_fn, opt); end



