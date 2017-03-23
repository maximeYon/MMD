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
    'delta','eta','s','p','l','fa','cs','cp','cl','cm'};
udtiparam = {'udelta','ucs','ucp','ucl','ufa'};
param = {dtiparam{:},udtiparam{:},'miso','viso','maniso','vaniso','msaniso','vsaniso'};
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
                    
                    udtd = [par'; perp'; 0*theta'; 0*phi'; w'];
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
dps.mu2macro = 4/45*((dps.lambda33-dps.lambda11).^2 + (dps.lambda22-dps.lambda11).*(dps.lambda22-dps.lambda33));

dps.mdelta = dps.maniso_n/3;
% dps.ufa = sqrt(3/2)*sqrt(1./((dps.miso.^2+dps.mu2iso)./(5/2*dps.mu2aniso)+1));
% dps.ufa(isnan(dps.ufa)) = 0;
% dps.ufa(dps.ufa>1) = 1;

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

dps.op = sqrt(dps.mu2macro./dps.mu2aniso);
dps.op(dps.op>1) = 1;
dps.op(isnan(dps.op)) = 0;

if (~isempty(dps_fn)) mdm_dps_save(dps, dps.s, dps_fn, opt); end



