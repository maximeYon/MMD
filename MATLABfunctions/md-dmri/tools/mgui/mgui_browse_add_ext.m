function [EG, c_ext] = mgui_browse_add_ext(EG, ext)
% function [EG, c_ext] = mgui_browse_add_ext(EG, ext)

if (any(isinf(ext)))
    ext = 'all';
    new_ext.id = ext;
    new_ext.ext = '';
    new_ext.name = 'All files';
    new_ext.filter = '*.*';
else
    new_ext.id = ext;
    new_ext.name = ext;
    new_ext.ext = ['.' ext];
    new_ext.filter = ['*.' ext];
end

if (~isfield(EG.browse,'ext'))
    EG.browse.ext{1} = new_ext;
end

% Don't add duplicates
for c_ext = 1:numel(EG.browse.ext)
    if (strcmp(EG.browse.ext{c_ext}.id, ext)), return; end
end

EG.browse.ext{end+1} = new_ext;
