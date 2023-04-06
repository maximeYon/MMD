function sz = msf_size(m, d)
% function sz = msf_size(m, d)
% 
% a size function that makes sure to return data with the correct
% dimensionality
%
% m - matrix
% d - dimensions, i.e. length of sz

sz = ones(1,d);
sz(1:numel(size(m))) = size(m);
