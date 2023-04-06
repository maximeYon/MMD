function chisqn = msf_chisqn(signal,signal_fit,ind)
% function chisqn = msf_chisqn(signal,signal_fit,ind)
%
% Calculate normalized chi-square from measured and fitted signals

if (nargin < 3)
    ind = 1:numel(signal);
end

nind = numel(ind);
chisqn = sum((signal_fit(ind) - signal(ind)).^2)/nind;