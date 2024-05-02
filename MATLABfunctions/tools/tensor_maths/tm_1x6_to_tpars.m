function t = tm_1x6_to_tpars(t_1x6)
% function t = tm_1x6_to_tpars(t_1x6)

if (size(t_1x6,2) ~= 6), error('wrong dimension of tensor'); end

tmp = cell(1, size(t_1x6,1));
for c = 1:size(t_1x6,1)
    tmp{c} = tm_3x3_to_tpars(tm_1x6_to_3x3(t_1x6(c,:)));
end

f = fieldnames(tmp{1});

t.present = 1;
for c = 1:numel(f)
    for k = 1:numel(tmp)
    t.(f{c})(k,:) = tmp{k}.(f{c});
    end
end
t = rmfield(t, 'present');



