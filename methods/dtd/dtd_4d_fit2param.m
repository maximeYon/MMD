function dps = dtd_4d_fit2param(mfs_fn, dps_fn, opt)
% function dps = dtd_4d_fit2param(mfs_fn, dps_fn, opt)

if (nargin < 2), dps_fn = []; end
if (nargin < 3), opt = []; end

opt = mdm_opt(opt);
dps = mdm_mfs_load(mfs_fn);

% create parameter maps and save them

sz = size(dps.m);
dps.s0 = zeros([sz(1) sz(2) sz(3)]);
dps.t1x6 = zeros([sz(1) sz(2) sz(3) 6]);
dps.lambdazzvec = zeros([sz(1) sz(2) sz(3) 3]);
dps.lambdaxxvec = zeros([sz(1) sz(2) sz(3) 3]);
dps.lambdayyvec = zeros([sz(1) sz(2) sz(3) 3]);
dps.lambda11vec = zeros([sz(1) sz(2) sz(3) 3]);
dps.lambda22vec = zeros([sz(1) sz(2) sz(3) 3]);
dps.lambda33vec = zeros([sz(1) sz(2) sz(3) 3]);
dtiparam = {'trace','iso','lambda33','lambda22','lambda11','lambdazz','lambdaxx','lambdayy','vlambda',...
    'delta','eta','s','p','l','fa','cs','cl','cp','cm'};
param = {dtiparam{:},'miso','viso','maniso','vaniso','mvlambda','vvlambda'};
for nparam = 1:numel(param)
    eval(['dps.' param{nparam} ' = zeros([sz(1) sz(2) sz(3)]);']);
end
for nk = 1:sz(3)
    for nj = 1:sz(2)
        for ni = 1:sz(1)
            %ni = 11; nj = 9; nk = 1;
            if dps.mask(ni,nj,nk)
                m = squeeze(dps.m(ni,nj,nk,:))';
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

                    dps.s0(ni,nj,nk) = s0;
                    dps.miso(ni,nj,nk) = miso;
                    dps.viso(ni,nj,nk) = viso;
                    dps.maniso(ni,nj,nk) = maniso;
                    dps.vaniso(ni,nj,nk) = vaniso;
                    dps.mvlambda(ni,nj,nk) = mvlambda;
                    dps.vvlambda(ni,nj,nk) = vvlambda;

                    [dtd_nx6,w] = dtd_dist2nx6w(dtd);
                    dt1x6 = (dtd_nx6'*w)'/s0;
                    dt3x3 = tm_1x6_to_3x3(dt1x6);

                    dt = tm_3x3_to_tpars(dt3x3);

                    dps.t1x6(ni,nj,nk,:) = dt.t1x6;
                    dps.lambdazzvec(ni,nj,nk,:) = dt.lambdazzvec;
                    dps.lambdaxxvec(ni,nj,nk,:) = dt.lambdaxxvec;
                    dps.lambdayyvec(ni,nj,nk,:) = dt.lambdayyvec;
                    dps.lambda11vec(ni,nj,nk,:) = dt.lambda11vec;
                    dps.lambda22vec(ni,nj,nk,:) = dt.lambda22vec;
                    dps.lambda33vec(ni,nj,nk,:) = dt.lambda33vec;
                    for nparam = 1:numel(dtiparam)
                        eval(['dps.' dtiparam{nparam} '(ni,nj,nk) = dt.' dtiparam{nparam} ';']);
                    end
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

kronecker = permute(repmat([1 1 1 0 0 0]',[1 size(dps.s0,3) size(dps.s0,2) size(dps.s0,1)]),[4 3 2 1]);
dps.s1x6 = (dps.t1x6./repmat(dps.miso,[1 1 1 6]) - kronecker)./repmat(dps.mdelta,[1 1 1 6])/2;
dps.s1x6prim = (2*dps.s1x6 + kronecker)/3;
dps.s1x6(isnan(dps.s1x6)) = 0;
dps.s1x6prim(isnan(dps.s1x6prim)) = 0;

dps.slambdaxx = (dps.lambdaxx./dps.miso - 1)./dps.mdelta/2;
dps.slambdaxx(isnan(dps.slambdaxx)) = 0;
dps.slambdayy = (dps.lambdayy./dps.miso - 1)./dps.mdelta/2;
dps.slambdayy(isnan(dps.slambdayy)) = 0;
dps.slambdazz = (dps.lambdazz./dps.miso - 1)./dps.mdelta/2;
dps.slambdazz(isnan(dps.slambdazz)) = 0;

dps.slambdaxxprim = (2*dps.slambdaxx + 1)/3;
dps.slambdayyprim = (2*dps.slambdayy + 1)/3;
dps.slambdazzprim = (2*dps.slambdazz + 1)/3;

dps.slambdas = zeros([sz(1) sz(2) sz(3) 3]);
dps.slambdas(:,:,:,1) = dps.slambdaxxprim;
dps.slambdas(:,:,:,2) = dps.slambdayyprim;
dps.slambdas(:,:,:,3) = dps.slambdazzprim;
dps.slambda11prim = min(dps.slambdas,[],4);
dps.slambda33prim = max(dps.slambdas,[],4);
dps = rmfield(dps,'slambdas');

dps.cc = dps.cm./dps.cmu;
dps.cc(isnan(dps.cc)) = 0;
dps.cc(dps.cc>1) = 1;
dps.op = sqrt(dps.cc);

if (~isempty(dps_fn)) mdm_dps_save(dps, dps.s, dps_fn, opt); end



