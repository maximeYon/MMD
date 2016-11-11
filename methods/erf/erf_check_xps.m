function erf_check_xps( xps )

f = {'b_delta', 'b', 'n'};

for c = 1:numel(f)
    assert( isfield( xps, f{c} ), ['xps.' f{c} ' required']);
end