function s_protocol = mdm_philips_parse_txt_protocol(prot_fn)
% function s_protocol = mdm_philips_parse_txt_protocol(prot_fn)

1;

% Read the file contents
if (~exist(prot_fn,'file')), error('file not found (%s)', prot_fn); end
fid = fopen(prot_fn);
str = fread(fid, inf, 'char'); 
fclose(fid);

% Remove all spaces after ';'
str2 = zeros(size(str));
on = 1; c2 = 1;
for c = 1:numel(str)
    if (on == 0) && (str(c) ~= ' '), on = 1; end
    if (on), str2(c2) = str(c); c2 = c2 + 1; end
    if (str(c) == ';'), on = 0; end
end
str = str2(1:(c2-1));

% Remove newlines unless preceeded by a ';'
str2 = zeros(size(str));
on = 1; c2 = 1;
for c = 1:numel(str)
    if (c > 2)
        if (str(c) == 13) && (str(c-1) ~= ';'), continue; end
        if (str(c) == 10) && (str(c-1) == 13) && (str(c-2) ~= ';'), continue; end
    end
    str2(c2) = str(c); c2 = c2 + 1;
end
str = str2(1:(c2-1));

m = strsplit(char(str'), char(10));

for c = 1:numel(m)
    
    if (numel(m) == 0), continue; end
    
    % Format the field name
    field_name = m{c}(1:(find(m{c} == '=', 1, 'first')-1));
    if (numel(field_name) == 0), continue; end
    
    field_name(1:(find(field_name ~= ' ', 1, 'first')-1)) = '_';
    field_name = strtrim(field_name);
    field_name = strrep(field_name, '-', '_');
    field_name = strrep(field_name, ' ', '_');
    field_name = strrep(field_name, '(', '');
    field_name = strrep(field_name, ')', '');
    field_name = strrep(field_name, '.', '');
    field_name = strrep(field_name, '/', '');
    field_name = strrep(field_name, '?', '');
    field_name = strrep(field_name, '>', '');
    field_name = strrep(field_name, '<', '');
    field_name = strrep(field_name, '%', '');
    field_name = strtrim(field_name);
    
    % Format the field value
    field_value = m{c}( (find(m{c} == '=', 1, 'first')+1):(find(m{c} == ';', 1, 'first')-1));
    
    field_value = strtrim(field_value);
    
    if (field_value(1) == '"') && (field_value(end) == '"')
        field_value = field_value(2:end-1);
    else
        m_value = regexp(field_value, '([^,]+)(?:[,]*)(?:\s*)', 'tokens');
        tmp = [];
        c_tmp = 1;
        for c_value = 1:numel(m_value)
            if (sum(m_value{c_value}{1} == '(') > 0);
                
                s = m_value{c_value}{1};
                n = str2double(s(...
                    (strfind(s, '(')+1):(strfind(s, ')')-1)));
                
                v = str2double(s((strfind(s,')')+1):end));
                
                tmp( c_tmp - 1 + (1:n)) = v;
                c_tmp = c_tmp + n;
            else
                tmp(c_tmp) = str2double(m_value{c_value}{1});
                c_tmp = c_tmp + 1;
            end
        end
        field_value = tmp;
    end
    
    try
        s_protocol.(['f' field_name]) = field_value;
    catch me
        disp(me.message);
    end
    
end
