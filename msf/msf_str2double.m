function o = msf_str2double(str, delim)
% function o = msf_str2double(str, delim)
%
% Splits 'str' according to 'delim' and converts to double

if (nargin < 2), delim = ' '; end

o = zeros(1, sum(str == delim) + 1);

for c = 1:numel(o)
    
    pos = strfind(str, delim);
    
    if any(numel(o) == [1 c])
        
        o(c) = str2double(str);
        
    else
        
        o(c) = str2double(str(1:pos(1)-1));
        
        str = str(pos(1)+1:end);
        
    end
end

