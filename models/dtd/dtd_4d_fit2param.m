function res = dtd_4d_fit2param(mfs_fn, o_fn, opt)
% function fn = dtd_4d_fit2param(mfs_fn, o_path, opt)

if (nargin < 3), opt = []; end

res = -1;
    
opt = mdm_opt(opt);
mfs = mdm_mfs_load(mfs_fn);
h   = mfs.nii_h; 

% create parameter maps and save them

sz = size(mfs.m);
mfs.s0 = zeros([sz(1) sz(2) sz(3)]);
mfs.t1x6 = zeros([sz(1) sz(2) sz(3) 6]);
mfs.lambdazzvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.lambdaxxvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.lambdayyvec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.lambda11vec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.lambda22vec = zeros([sz(1) sz(2) sz(3) 3]);
mfs.lambda33vec = zeros([sz(1) sz(2) sz(3) 3]);
dtiparam = {'trace','iso','lambda33','lambda22','lambda11','lambdazz','lambdaxx','lambdayy','vlambda',...
    'delta','eta','s','p','l','fa','cs','cl','cp','cm'};
param = {dtiparam{:},'miso','viso','maniso','vaniso','mvlambda','vvlambda'};
for nparam = 1:numel(param)
    eval(['mfs.' param{nparam} ' = zeros([sz(1) sz(2) sz(3)]);']);
end
for nk = 1:sz(3)
    for nj = 1:sz(2)
        for ni = 1:sz(1)
            %ni = 11; nj = 9; nk = 1;
            if mfs.mask(ni,nj,nk)
                m = squeeze(mfs.m(ni,nj,nk,:))';
                dtd = dtd_m2dtd(m);
                [n,par,perp,theta,phi,w] = dtd_dist2par(dtd);
                                
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

                    [dtd_nx6,w] = dtd_dist2nx6w(dtd);
                    dt1x6 = (dtd_nx6'*w)'/s0;
                    dt3x3 = tm_1x6_to_3x3(dt1x6);

                    dt = tm_t2tpars(dt3x3);

                    mfs.t1x6(ni,nj,nk,:) = dt.t1x6;
                    mfs.lambdazzvec(ni,nj,nk,:) = dt.lambdazzvec;
                    mfs.lambdaxxvec(ni,nj,nk,:) = dt.lambdaxxvec;
                    mfs.lambdayyvec(ni,nj,nk,:) = dt.lambdayyvec;
                    mfs.lambda11vec(ni,nj,nk,:) = dt.lambda11vec;
                    mfs.lambda22vec(ni,nj,nk,:) = dt.lambda22vec;
                    mfs.lambda33vec(ni,nj,nk,:) = dt.lambda33vec;
                    for nparam = 1:numel(dtiparam)
                        eval(['mfs.' dtiparam{nparam} '(ni,nj,nk) = dt.' dtiparam{nparam} ';']);
                    end
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

kronecker = permute(repmat([1 1 1 0 0 0]',[1 size(mfs.s0,3) size(mfs.s0,2) size(mfs.s0,1)]),[4 3 2 1]);
mfs.s1x6 = (mfs.t1x6./repmat(mfs.miso,[1 1 1 6]) - kronecker)./repmat(mfs.mdelta,[1 1 1 6])/2;
mfs.s1x6prim = (2*mfs.s1x6 + kronecker)/3;
mfs.s1x6(isnan(mfs.s1x6)) = 0;
mfs.s1x6prim(isnan(mfs.s1x6prim)) = 0;

mfs.slambdaxx = (mfs.lambdaxx./mfs.miso - 1)./mfs.mdelta/2;
mfs.slambdaxx(isnan(mfs.slambdaxx)) = 0;
mfs.slambdayy = (mfs.lambdayy./mfs.miso - 1)./mfs.mdelta/2;
mfs.slambdayy(isnan(mfs.slambdayy)) = 0;
mfs.slambdazz = (mfs.lambdazz./mfs.miso - 1)./mfs.mdelta/2;
mfs.slambdazz(isnan(mfs.slambdazz)) = 0;

mfs.slambdaxxprim = (2*mfs.slambdaxx + 1)/3;
mfs.slambdayyprim = (2*mfs.slambdayy + 1)/3;
mfs.slambdazzprim = (2*mfs.slambdazz + 1)/3;

mfs.slambdas = zeros([sz(1) sz(2) sz(3) 3]);
mfs.slambdas(:,:,:,1) = mfs.slambdaxxprim;
mfs.slambdas(:,:,:,2) = mfs.slambdayyprim;
mfs.slambdas(:,:,:,3) = mfs.slambdazzprim;
mfs.slambda11prim = min(mfs.slambdas,[],4);
mfs.slambda33prim = max(mfs.slambdas,[],4);
mfs = rmfield(mfs,'slambdas');

mfs.cc = mfs.cm./mfs.cmu;
mfs.cc(isnan(mfs.cc)) = 0;
mfs.cc(mfs.cc>1) = 1;
mfs.op = sqrt(mfs.cc);

mfs_fn = mdm_mfs_save(mfs, mfs.s, o_fn, opt);

res = 1;

