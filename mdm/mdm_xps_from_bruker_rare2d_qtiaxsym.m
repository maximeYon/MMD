function res = mdm_xps_from_bruker_rare2d_qtiaxsym(data_path)
% function res = mdm_xps_from_bruker_rare2d_dti(data_path)
%
% read gradient in the directory data_path
% save xps
%
% the xps is always in SI units

load([data_path '/NMRacqus'])

if any(strcmp(NMRacqus.pulprog,{'DT_qVASrare2d'})) == 1

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

    load([data_path '/DiffRamp'])

    NbTtrace = DiffRamp.NbTtrace;
    NbTDelta = DiffRamp.NbTDelta;
    NbTdir = DiffRamp.NbTdir;
    NbTtot = NbTtrace*NbTDelta*NbTdir;
    GridSamp = DiffRamp.GridSamp;
    td1 = NbTtot;

    fid = fopen([data_path '/rax.txt']);
    ramp.ax = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/ray.txt']);
    ramp.ay = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/raz.txt']);
    ramp.az = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/rbx.txt']);
    ramp.bx = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/rby.txt']);
    ramp.by = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/rbz.txt']);
    ramp.bz = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/rcx.txt']);
    ramp.cx = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/rcy.txt']);
    ramp.cy = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/rcz.txt']);
    ramp.cz = fscanf(fid,'%f');
    fclose(fid);

    G.ax = ramp.ax*Gmax.*NMRacqus.cnst1/100;
    G.bx = ramp.bx*Gmax.*NMRacqus.cnst1/100;
    G.cx = ramp.cx*Gmax.*NMRacqus.cnst1/100;
    G.ay = ramp.ay*Gmax.*NMRacqus.cnst2/100;
    G.by = ramp.by*Gmax.*NMRacqus.cnst2/100;
    G.cy = ramp.cy*Gmax.*NMRacqus.cnst2/100;
    G.az = ramp.az*Gmax.*NMRacqus.cnst3/100;
    G.bz = ramp.bz*Gmax.*NMRacqus.cnst3/100;
    G.cz = ramp.cz*Gmax.*NMRacqus.cnst3/100;

    %symmetry vector of b-matrix
    symv.x = G.cx;
    symv.y = G.cy;
    symv.z = G.cz;
    symv.norm = sqrt(symv.x.^2 + symv.y.^2 + symv.z.^2);
    symv.x = symv.x./symv.norm;
    symv.y = symv.y./symv.norm;
    symv.z = symv.z./symv.norm;
    symv.theta = acos(symv.z);
    symv.phi = atan2(symv.y,symv.x);

    %qVAS gradient time-modulation
    fid = fopen([data_path '/qVASa.txt']);
    Gmod.a = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/qVASb.txt']);
    Gmod.b = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/qVASc.txt']);
    Gmod.c = fscanf(fid,'%f');
    fclose(fid);

    Ndt = length(Gmod.c);
    tau = NMRacqus.d3;
    dt = tau/Ndt;

    %figure(1), clf, plot((1:Ndt)',[Gmod.c Gmod.b Gmod.a],'-'), return

    G.x = repmat(Gmod.a,[1, td1]).*repmat(G.ax',[Ndt, 1]) + ...
        repmat(Gmod.b,[1, td1]).*repmat(G.bx',[Ndt, 1]) + ...
        repmat(Gmod.c,[1, td1]).*repmat(G.cx',[Ndt, 1]);
    G.y = repmat(Gmod.a,[1, td1]).*repmat(G.ay',[Ndt, 1]) + ...
        repmat(Gmod.b,[1, td1]).*repmat(G.by',[Ndt, 1]) + ...
        repmat(Gmod.c,[1, td1]).*repmat(G.cy',[Ndt, 1]);
    G.z = repmat(Gmod.a,[1, td1]).*repmat(G.az',[Ndt, 1]) + ...
        repmat(Gmod.b,[1, td1]).*repmat(G.bz',[Ndt, 1]) + ...
        repmat(Gmod.c,[1, td1]).*repmat(G.cz',[Ndt, 1]);

    %dephasing vector F
    F.x = cumsum(G.x*dt,1);
    F.y = cumsum(G.y*dt,1);
    F.z = cumsum(G.z*dt,1);
    F.r = sqrt(F.x.^2 + F.y.^2 + F.z.^2);

    %diffusion weighting matrix b
    %factor 2 from the double DW blocks
    bt.xx = 2*gamma^2*sum(F.x.*F.x*dt,1)';
    bt.xy = 2*gamma^2*sum(F.x.*F.y*dt,1)';
    bt.xz = 2*gamma^2*sum(F.x.*F.z*dt,1)';
    bt.yy = 2*gamma^2*sum(F.y.*F.y*dt,1)';
    bt.yz = 2*gamma^2*sum(F.y.*F.z*dt,1)';
    bt.zz = 2*gamma^2*sum(F.z.*F.z*dt,1)';

    b = (bt.xx + bt.yy + bt.zz);
    %figure(1), clf, plot((1:td1)',[b],'o'), return
    
    xps.b = b;
    xps.n = td1;
    xps.bt = [bt.xx bt.yy bt.zz sqrt(2)*[bt.xy bt.xz bt.yz]];
    xps.gma = Gmod.a';
    xps.gmb = Gmod.b';
    xps.gmc = Gmod.c';
    xps.tau = tau;
    xps.ngm = 2;
    
    xps.b_nind = NbTtrace;
    xps.bd_nind = NbTDelta;
    xps.br_nind = NbTdir;
    [b_ind,bd_ind,br_ind] = ndgrid(1:xps.b_nind,1:xps.bd_nind,1:xps.br_nind);
    xps.b_ind = reshape(b_ind,td1,1);
    xps.bd_ind = reshape(bd_ind,td1,1);
    xps.br_ind = reshape(br_ind,td1,1);

    xps_old = xps;
    xps = mdm_xps_btpars_from_bt(xps_old);

end

save([data_path '/xps'],'xps')

res = 1;
