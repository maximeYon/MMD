function res = dtd_pa_4d_fit2param(mfs_fn, o_fn, opt)
% function fn = dtd_4d_fit2param(mfs_fn, o_path, opt)

if (nargin < 3), opt = []; end

res = -1;
    
opt = mdm_opt(opt);
mfs = mdm_mfs_load(mfs_fn);
h   = mfs.nii_h; 

% create parameter maps and save them

sz = size(mfs.m);
param = {'s0','miso','viso','maniso','vaniso','mvlambda','vvlambda'};
for nparam = 1:numel(param)
    eval(['mfs.' param{nparam} ' = zeros([sz(1) sz(2) sz(3)]);']);
end
for nk = 1:sz(3)
    for nj = 1:sz(2)
        for ni = 1:sz(1)
            %ni = 11; nj = 9; nk = 1;
            if mfs.mask(ni,nj,nk)
                m = squeeze(mfs.m(ni,nj,nk,:))';
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

                    mfs.s0(ni,nj,nk) = s0;
                    mfs.miso(ni,nj,nk) = miso;
                    mfs.viso(ni,nj,nk) = viso;
                    mfs.maniso(ni,nj,nk) = maniso;
                    mfs.vaniso(ni,nj,nk) = vaniso;
                    mfs.mvlambda(ni,nj,nk) = mvlambda;
                    mfs.vvlambda(ni,nj,nk) = vvlambda;
                end                
            end
        end
    end
end

mfs.mu2iso = mfs.viso;
mfs.mu2aniso = 2/5*mfs.mvlambda;
mfs.mu2aniso(isnan(mfs.mu2aniso)) = 0;
mfs.mu2tot = mfs.mu2iso + mfs.mu2aniso;

mfs.kiso = mfs.mu2iso./mfs.miso.^2;
mfs.kiso(isnan(mfs.kiso)) = 0;
mfs.kaniso = mfs.mu2aniso./mfs.miso.^2;
mfs.kaniso(isnan(mfs.kaniso)) = 0;
mfs.ktot = mfs.mu2tot./mfs.miso.^2;
mfs.ktot(isnan(mfs.ktot)) = 0;

mfs.ciso = mfs.viso./mfs.miso.^2;
mfs.ciso(isnan(mfs.ciso)) = 0;
mfs.mdelta = mfs.maniso./(3*mfs.miso);
mfs.mdelta(isnan(mfs.mdelta)) = 0;
mfs.vdelta = mfs.vaniso./(3*mfs.miso).^2;
mfs.vdelta(isnan(mfs.vdelta)) = 0;

mfs.cmvlambda = mfs.mvlambda./mfs.miso.^2;
mfs.cmvlambda(isnan(mfs.cmvlambda)) = 0;
mfs.ufa = sqrt(3/2)*sqrt(1./(1./mfs.cmvlambda+1));
mfs.ufa(isnan(mfs.ufa)) = 0;
mfs.cmu = mfs.ufa.^2;

mfs.cvvlambda = mfs.vvlambda./mfs.miso.^4;
mfs.cvvlambda(isnan(mfs.cvvlambda)) = 0;

mfs_fn = mdm_mfs_save(mfs, mfs.s, o_fn, opt);

res = 1;

