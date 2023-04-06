function chisq_fn = mio_datafit2chisq(fun, s, mfs_fn, o_fn, opt)
% function chisq_fn = mio_datafit2chisq(fun, s, mfs_fn, o_fn, opt)
%

% Read and reformat data
[I,h]  = mdm_nii_read(s.nii_fn);
M      = mdm_mask_load(s, opt);
xps = s.xps;
mfs = mdm_mfs_load(mfs_fn);
m = mfs.m;
sz = size(m);

ind = opt.dtd.ind_start:s.xps.n;
nind = numel(ind);

% Disallow model fits to complex data
if (any(imag(I) ~= 0)), I = abs(I); end 

chisq = zeros([sz(1) sz(2) sz(3)]);

% Loop over the data
for nk = 1:sz(3)
    for nj = 1:sz(2)        
        for ni = 1:sz(1)
            if M(ni,nj,nk)
                signal_fit = fun(squeeze(m(ni,nj,nk,:)), xps);
                signal = double(squeeze(I(ni,nj,nk,:)));
                chisq(ni,nj,nk) = sum((signal_fit(ind) - signal(ind)).^2)/nind;
            end
        end
    end
end

save(o_fn, 'chisq')
chisq_fn = o_fn;