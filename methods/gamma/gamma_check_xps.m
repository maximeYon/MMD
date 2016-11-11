function gamma_check_xps( xps )
% function gamma_check_xps( xps )

f = {'b_delta', 'b_eta', 'b', 'n'};

for c = 1:numel(f)
    assert( isfield( xps, f{c} ), ['xps.' f{c} ' required']);
end