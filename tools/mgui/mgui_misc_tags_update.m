function tags = mgui_misc_tags_update(tags, new_tags)
% function tags = mgui_misc_tags_update(tags, new_tags)

for c_tag = 1:numel(new_tags)
    try
        if (~isfield(tags, new_tags{c_tag}) || isempty(tags.(new_tags{c_tag})))
            tags.(new_tags{c_tag}) = findobj('Tag', new_tags{c_tag});
        end
    catch me
        disp(me.message);
    end
end
