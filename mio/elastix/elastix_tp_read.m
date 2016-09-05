function [p,tmat,str] = elastix_tp_read(fn)
% function [p,tmat,str] = elastix_tp_read(fn)
%
% Reads a transform parameter file

% Read the whole file
fid = fopen(fn);
str = fread(fid, inf, 'uint8=>char')';
fclose(fid);

if isempty(strfind(lower(str), lower('(Transform "AffineDTITransform")')))
    error('elastix_tp_read: for now, this function only works for AffineDTITransform');
end

% Find the TransformParameters line
t = regexp(str, '(?:TransformParameters[ ]+)([^)]*)', 'tokens');

if (numel(t) ~= 1), error('Multiple transformparameters in one file?'); end

t = regexp(t{1}, '([^ ]*)(?:[ ]*)', 'tokens');
t = t{1};

p = zeros(1,size(t,2));
for c = 1:size(t,2)
    p(c) = str2double(t{c});
end

% Reorder and reshape the parameters as rot, tra, scale, skew (as Leemans)
p = reshape(p([1 2 3   10 11 12    7 8 9    4 5 6]), 3, 4);

% Also read the translation matrix
if (nargout > 1)
    t = regexp(str, '(?:MatrixTranslation[ ]+)([^)]*)', 'tokens');
    
    if (numel(t) ~= 1), error('Multiple transformparameters in one file?'); end
    
    t = regexp(t{1}, '([^ ]*)(?:[ ]*)', 'tokens');
    t = t{1};
    
    tmat = zeros(1,size(t,2));
    for c = 1:size(t,2)
        tmat(c) = str2double(t{c});
    end
    tmat = [reshape(tmat, [3 4]); 0 0 0 1];
    
end


