function gamma = msf_const_gamma(nucleus)
% function gamma = msf_const_gamma(nucleus)

if (nargin < 1), nucleus = '1H'; end

if (~strcmp(nucleus, '1H')), error('not implemented'); end

gamma = 42.576e6 * 2 * pi;
