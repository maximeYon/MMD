function m = dtd_dtd2m(dtd,opt)

m = zeros(1 + 5*opt.dtd.n_out,1);
if ~isempty(dtd)
    m(1:numel(dtd)) = dtd;
end
m = m';
