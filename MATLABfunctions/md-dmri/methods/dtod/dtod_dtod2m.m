function m = dtod_dtod2m(dtod,opt)

m = zeros(1 + 8*opt.dtod.n_out,1);
if ~isempty(dtod)
    m(1:numel(dtod)) = dtod;
    m = m(1:(1 + 8*opt.dtod.n_out),1);
    m(1) = min([dtod(1) opt.dtod.n_out]);
end
m = m';
