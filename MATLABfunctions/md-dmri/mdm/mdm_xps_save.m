function mdm_xps_save(xps, xps_fn, opt)
% function mdm_xps_save(xps, xps_fn)

if (nargin < 3), opt.present = 1; end
opt = mdm_opt(opt);

if (exist(xps_fn, 'file') && (~opt.do_overwrite))
    msf_fprintf(opt, 'File found, skipping (%s)\n', xps_fn);
    return;
end

save(xps_fn, 'xps');


