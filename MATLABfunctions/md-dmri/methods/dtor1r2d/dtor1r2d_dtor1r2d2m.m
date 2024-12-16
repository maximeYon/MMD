function m = dtor1r2d_dtor1r2d2m(dtor1r2d,opt)

m = zeros(1 + 10*opt.dtor1r2d.n_out,1);
if ~isempty(dtor1r2d)
    m(1:numel(dtor1r2d)) = dtor1r2d;
    m = m(1:(1 + 10*opt.dtor1r2d.n_out),1);
    m(1) = min([dtor1r2d(1) opt.dtor1r2d.n_out]);
end
m = m';
