function [A, xps_pa] = mio_powder_average(I, xps, opt)
% function [J, xps] = mio_powder_average(I, xps)

if (nargin < 3), opt = []; end

opt = mdm_opt(opt);

I = double(I);

% Average image
id = xps.a_ind;
if (isfield(xps,'s_ind'))
    id = [id xps.s_ind];
end
[~,~,id_ind] = unique(id, 'rows');

% get rid of NaNs
tmp = sum(isnan(id),2) > 0;
c_list = unique(id_ind(~tmp)); 
n = numel(c_list);

A = zeros(size(I,1), size(I,2), size(I,3), n);

for c = c_list'
    A(:,:,:,c == c_list) = msf_nanmean(I(:,:,:,id_ind == c),4);
end

if (opt.do_pa_abs), A = abs(A); end

% Average fields in xps
xps = rmfield(xps, 'n');
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
            if (opt.pa_rethrow_error)
                error('failed powder averaging field %s', f{i});
            else
                warning('failed powder averaging field %s', f{i});
            end
        end
    end
end
xps_pa.n = n;

