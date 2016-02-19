function res = mdm_xps_from_bruker_rare2d_dti(data_path)
% function res = mdm_xps_from_bruker_rare2d_dti(data_path)
%
% read gradient in the directory data_path
% save xps
%
% the xps is always in SI units

eval(['load ' data_path '/NMRacqus'])

if any(strcmp(NMRacqus.pulprog,{'DT_dtirare2d','DT_dpgsedtirare2d'})) == 1

    gamma = 26.75e7;
    if strcmp(NMRacqus.nuc1,'2H') == 1
        gamma = 4.1065e7;
    elseif strcmp(NMRacqus.nuc1,'23Na') == 1
        gamma = 7.0761e7;
    end

    Gmax = 3;
    if any(strcmp(NMRacqus.probhd,{'5 mm BBO BB-1H/D XYZ-GRD Z107255/0001',...
            '5 mm TXI 1H/D-13C/15N XYZ-GRD Z8588/0006'})) == 1
        Gmax = 0.5;
    end

    fid = fopen([data_path '/difframp_x.txt']);
    difframp.x = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/difframp_y.txt']);
    difframp.y = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/difframp_z.txt']);
    difframp.z = fscanf(fid,'%f');
    fclose(fid);

    G.x = difframp.x*Gmax.*NMRacqus.cnst1/100;
    G.y = difframp.y*Gmax.*NMRacqus.cnst2/100;
    G.z = difframp.z*Gmax.*NMRacqus.cnst3/100;
    G.r = sqrt(G.x.^2 + G.y.^2 + G.z.^2);
    Gnorm.x = G.x./G.r;
    Gnorm.y = G.y./G.r;
    Gnorm.z = G.z./G.r;

    if any(strcmp(NMRacqus.pulprog,{'DT_dtirare2d'})) == 1
        delta = NMRacqus.d3+NMRacqus.d2; epsilon = NMRacqus.d2;
        Delta = delta+epsilon+2*NMRacqus.d4+NMRacqus.p12/1e6+4*NMRacqus.d32;
        tdiff = Delta-delta/3+epsilon^3/30/delta^2-epsilon^2/6/delta;
        q = gamma*G.r*delta/2/pi;
        b = (2*pi*q).^2*tdiff;
    elseif any(strcmp(NMRacqus.pulprog,{'DT_dpgsedtirare2d'})) == 1
        delta = NMRacqus.d3+NMRacqus.d2; epsilon = NMRacqus.d2;
        Delta = delta+epsilon+2*NMRacqus.d4+NMRacqus.p12/1e6+4*NMRacqus.d32;
        tdiff = Delta-delta/3+epsilon^3/30/delta^2-epsilon^2/6/delta;
        q = gamma*G.r*delta/2/pi;
        b = 2*(2*pi*q).^2*tdiff;
    end

    xps.b   = b; 
    xps.u   = [Gnorm.x Gnorm.y Gnorm.z];
    xps.n   = numel(xps.b);
    xps.bt  = dtd_1x3_to_1x6(xps.b, zeros(size(xps.b)), xps.u);
    xps.bt2 = dtd_1x6_to_1x21(xps.bt);
end

save([data_path '/xps'],'xps')

res = 1;
