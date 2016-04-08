function xps = mdm_xps_from_bval_bvec(bval_fn, bvec_fn)
% function xps = mdm_xps_from_bval_bvec(bval_fn, bvec_fn)

bval = mdm_txt_read(bval_fn);
xps.b = str2num(bval{1})' * 1e6; 
xps.n = numel(xps.b);


bvec = mdm_txt_read(bvec_fn);
assert(numel(bvec) == 3, 'strange bvec file');


gdir = [str2num(bvec{1}); str2num(bvec{2}); str2num(bvec{3})];

xps.u = gdir';

xps.bt = tm_1x3_to_1x6(xps.b, zeros(size(xps.b)), xps.u);
