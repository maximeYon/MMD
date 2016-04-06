function m = dtd_pa_dtd2m(dtd,opt)

m = zeros(1 + 3*opt.dtd_pa.n_out,1);
if ~isempty(dtd)
    m(1:numel(dtd)) = dtd;
end
m = m';
