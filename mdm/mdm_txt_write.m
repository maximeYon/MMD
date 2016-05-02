function txt = mdm_txt_write(txt, txt_fn, opt)
% function txt = mdm_txt_write(txt, txt_fn, opt)
% 
% This function write each line in txt, a cell, as separate lines in txt_fn

if (nargin < 3), opt.present = 1; end

opt = mdm_opt(opt);

if (~iscell(txt)), error('txt should be cell'); end

fid = fopen(txt_fn, 'w');
for c = 1:numel(txt)
    fprintf(fid, '%s\n', txt{c});
end
fclose(fid);
