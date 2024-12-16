function dps_fn = mdm_fit2param(fun_4d_fit2param, mfs_fn, dps_fn, opt)
% function dps_fn = mdm_fit2param(fun_4d_fit2param, mfs_fn, dps_fn, opt)
%
% fun_4d_fit2param - parameter derivation function
% mfs_fn           - model fit structure read from dps_fn
% dps_fn           - derived parameter structure is written to dps_fn
% opt              - options structure

if (nargin < 4), opt.present = 1; end

opt = mdm_opt(opt);

msf_log(['Starting ' mfilename], opt);    

if (~iscell(mfs_fn))
    assert(~strcmp(mfs_fn, dps_fn), 'pats.mfs_fn and dps_fn must be different');
else
    for c = 1:numel(mfs_fn)
        assert(~strcmp(mfs_fn{c}, dps_fn), 'pats.mfs_fn{%i} and dps_fn must be different', c);
    end
end

if (exist(dps_fn, 'file') && (~opt.do_overwrite))
    disp(['Skipping, output file already exists: ' dps_fn]); return;
end

if (~opt.do_fit2param)
    dps_fn = '';
    return;
end

dps_fn = fun_4d_fit2param(mfs_fn, dps_fn, opt);


