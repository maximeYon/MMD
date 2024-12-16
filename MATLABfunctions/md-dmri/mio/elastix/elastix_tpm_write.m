function elastix_tpm_write(tpm, fn)
% function elastix_tpm_write(tpm, fn)
%
% Writes the 'transform parameter matrix' to the file named 'fn'

fid = fopen(fn, 'w+');
for c = 1:size(tpm,2)
    fprintf(fid, [sprintf('%1.6f ', tpm(:,c)) '\n']);
end
fclose(fid);
