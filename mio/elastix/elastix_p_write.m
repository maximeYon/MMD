function elastix_param_write(p, fn)
%function elastix_param_write(p, fn)

if (~isstruct(p))
    error('expected first argument to be a structure');
end

f = fieldnames(p);

str = '';
for c = 1:numel(f)
    
    if (isnumeric(p.(f{c})))
        tmp = p.(f{c});
        if (numel(tmp) == 1) && ( round(round(tmp) - tmp) < eps)
            value = sprintf('%i ', tmp(:)');
        else
            value = sprintf('%f ', tmp(:)');
        end
    elseif (any(p.(f{c}) == '"'))
        value = p.(f{c});
    else
        value = ['"' p.(f{c})  '"'];
    end
    
    str = [str sprintf(['(' f{c} ' ' value ')\n'])]; %#ok<AGROW>
end

fid = fopen(fn, 'w');
fwrite(fid, char(str), 'char');
fclose(fid);

    