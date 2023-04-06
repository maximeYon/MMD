function dps = mio_1d_fit2param(method_name,m)
% function dps = mio_1d_fit2param(method_name,m)

mfs.m = zeros(1,1,1,numel(m)); mfs.m(1,1,1,:) = m;
dps = feval([method_name '_4d_fit2param'],mfs.m);

