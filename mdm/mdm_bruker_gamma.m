function gamma = mdm_bruker_gamma(NMRacqus)
% function gamma = mdm_bruker_gamma(NMRacqus)
%
% NMRacqus: Bruker acquistion parameter structure
% NMRacqus.nuc1: nucleus in acquisiton channel 1

if strcmp(NMRacqus.nuc1,'1H') == 1
   gamma = 26.75e7; %rad/(T*s)
elseif strcmp(NMRacqus.nuc1,'2H') == 1
    gamma = 4.1065e7;
elseif strcmp(NMRacqus.nuc1,'23Na') == 1
    gamma = 7.0761e7;
end
