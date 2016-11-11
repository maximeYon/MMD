function ndi16_check_xps( xps )
% function ndi16_check_xps( xps )
%
% required xps-fields: b and n

f = {'b', 'n'};

for c = 1:numel(f)
    assert( isfield( xps, f{c} ), ['xps.' f{c} ' required']);
end

