function [signal_pa,xps_pa] = mdm_powder_average_1d(signal,xps)
% function s = mdm_powder_average(s, o_path, opt)
%
% Average over rotations. Image volumes with identical rotations is defined
% from s.xps.a_ind
%
% To do: find a way of keeping track of number of averages per step

id = xps.a_ind;
[~,~,id_ind] = unique(id, 'rows');

% get rid of NaNs
tmp = sum(isnan(id),2) > 0;
c_list = unique(id_ind(~tmp)); 
n = numel(c_list);

signal_pa = zeros(n,1);

for c = c_list'
    signal_pa(c == c_list,1) = nanmean(signal(id_ind == c,1),1);
end

% Average fields in xps
xps_pa.n = n;
xps_pa.n_bs = xps.n_bs;
xps_pa.n_bl = xps.n_bl;

xps = rmfield(xps, {'n','n_bs','n_bl','n_br'});
f = fieldnames(xps);
for i = 1:numel(f)
    for c = c_list'

        % allow text fields to be just copied
        if (all(ischar(xps.(f{i}))))
            xps_pa.(f{i}) = xps.(f{i});
            continue; 
        end

        try
            xps_pa.(f{i})(c == c_list,:) = mean(xps.(f{i})(id_ind == c, :), 1);
        catch
            error('failed powder averaging field %s', f{i});
        end
    end
end
