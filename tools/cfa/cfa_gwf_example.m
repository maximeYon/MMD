function [gwf, rf, dt] = cfa_gwf_example(mode)
% function [gwf, rf, dt] = cfa_gwf_example(mode)
% Container for some example waveforms in the CFA framework.

switch mode
    case 1 % SE-SDE in x
        gwf = [
            0 1 1 1 1 1 0 0 0 1 1 1 1 1 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            ]'*.08;
        
        rf  = [1 1 1 1 1 1 0 0 0 -1 -1 -1 -1 -1 -1]';
        
        dt = 3e-3;
        
    case 2 % SE-DDE in x and y, respectively (no error if not rotated)
        gwf = [
            0 1 1 0 -1 -1 0 0 0 0 0 0  0  0 0;
            0 0 0 0  0  0 0 0 0 1 1 0 -1 -1 0;
            0 0 0 0  0  0 0 0 0 0 0 0  0  0 0
            ]'*.08;
        
        rf  = [1 1 1 1 1 1 0 0 0 -1 -1 -1 -1 -1 -1]';
        
        dt = 6e-3;
        
    case 3 % SE-DDE in x and z, respectively (significant error w/o rot)
        gwf = [
            0 1 1 0 -1 -1 0 0 0 0 0 0  0  0 0;
            0 0 0 0  0  0 0 0 0 0 0 0  0  0 0;
            0 0 0 0  0  0 0 0 0 1 1 0 -1 -1 0;
            ]'*.08;
        
        rf  = [1 1 1 1 1 1 0 0 0 -1 -1 -1 -1 -1 -1]';
        
        dt = 6e-3;
end

