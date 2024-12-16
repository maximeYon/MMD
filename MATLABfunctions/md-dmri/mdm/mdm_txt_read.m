function txt = mdm_txt_read(txt_fn, opt)
% function txt = mdm_txt_read(txt_fn)
% 
% This function returns each line in txt_fn as a cell, skipping empty
% lines and those starting with #
%

if (nargin < 2), opt.present = 1; end

opt = mdm_opt(opt);

fid = fopen(txt_fn);
txt = {};
while (1)
    
    % get line from file, break if end of file
    tline = fgetl(fid);
    if (~ischar(tline)), break, end
    tline = strtrim(tline);
    
    % skip the line under some circumstances
    if (numel(tline) == 0), continue; end
    if ((tline(1) == '#') && opt.mdm.txt_read_skip_comments), continue; end
    
    txt{end+1} = tline;
    
end
fclose(fid);