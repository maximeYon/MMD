function codivide16_check_xps( xps )
% function codivide16_check_xps( xps )
%
% xps-fields required: b_delta, b, and n

f = {'b_delta', 'b', 'n'};

for c = 1:numel(f)
    assert( isfield( xps, f{c} ), ['xps.' f{c} ' required']);
end

