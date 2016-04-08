function gamma = mdm_bruker_gamma(NMRacqus)
% function gamma = mdm_bruker_gamma(NMRacqus)
%
% NMRacqus: Bruker acquistion parameter structure
% NMRacqus.nuc1: nucleus in acquisiton channel 1

switch (strtrim(NMRacqus.nuc1))
    case '1H'
        gamma = 26.75e7; %rad/(T*s)
    case '2H'
        gamma = 4.1065e7;
    case '23Na'
        gamma = 7.0761e7;
end

