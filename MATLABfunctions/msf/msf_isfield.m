function val = msf_isfield(s, field_name)
% function val = msf_isfield(s, field_name)

if (verLessThan('matlab','8.4.0'))
    val = isfield(s, field_name);
else
    switch (class(s))
        case 'struct'
            val = isfield(s, field_name);
        otherwise
            val = isprop(s, field_name);
    end
end

