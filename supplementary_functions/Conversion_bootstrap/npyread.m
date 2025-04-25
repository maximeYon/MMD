function X = npyread(F)
% X = npyread(F)
%   Read an array from an npy file into a Matlab array.
%   Not universal, but hacked to handle the ones we need to cover.
%   See
%   https://github.com/numpy/numpy/blob/master/doc/neps/npy-format.txt
%
% 2013-05-10 Dan Ellis dpwe@ee.columbia.edu

if ~exist(F, 'file') 
  error(['npyread: no such file: ',F]);
end

fp = fopen(F, 'rb');

% Assume we're OK now
ver = fread(fp, 2, 'uchar');

header_len = fread(fp, 1, 'short');

fmt = char(fread(fp, header_len, 'char')');

% parse
t = regexp(fmt, '''(?<key>\w+)'' *: *''(?<val>[^'']*)''|''(?<key>\w+)'' *: *\((?<val>.*)\)|''(?<key>\w+)'' *: *(?<val>\w*)', 'names');

shape = [];

for i = 1:length(t)
  if strcmp(t(i).key, 'descr')
    if strcmp(t(i).val, '|O')
      format = 'float32';
    else
      fclose(fp);
      error(['npyread: unrecognized ''descr'' ''', t(i).val, [''' is  not ']'<f4''']);
    end
  elseif strcmp(t(i).key, 'fortran_order')
    if strcmp(t(i).val, 'True');
      fortran_order = 1;
    elseif strcmp(t(i).val, 'False');
      fortran_order = 0;
    else
      fclose(fp);
      error(['npyread: fortran_order = ',t(i).val,' is not True or False']);
    end
  elseif strcmp(t(i).key, 'shape');
    s = t(i).val;
    if s(end) == ','
      s = s(1:end-1);
    end
    shape = cellfun(@str2num, regexp(s, ' *, *', 'split'));
  else
    disp(['npyread: unrecognized format field ', t(i).key, ' : ', ...
          t(i).val]);
  end
end

% Assume shape is set and that it's a 2D array
if length(shape) == 1
  shape = [shape,1];
end

assert(length(shape) == 2);

assert(strcmp(format, 'float32'));
X = fread(fp, shape(2)*shape(1), 'float32');

X = reshape(X, shape(2), shape(1));

fclose(fp);