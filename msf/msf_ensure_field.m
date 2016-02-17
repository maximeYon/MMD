function s = msf_ensure_field(s, f, v)

if (~isfield(s, f)), s.(f) = v; end

