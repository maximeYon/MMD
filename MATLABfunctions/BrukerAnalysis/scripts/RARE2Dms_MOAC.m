clear allwd = cd;%DataDir = '/opt/topspin2/data/DT/nmr';%DataDir = '/Users/daniel/NMRdata/AVII500/DT';%DataDir = '/Users/daniel/Documents/Spaces/Presentations';%DataDir = '/Users/daniel/Dropbox';DataDir = '/Users/daniel/Dropbox/NMRdata/AVII500';ExpNam = {'MOACRARE2D_test'}; expno = 90;%lb = 150e-6; si = 128; si1 = si;%lb = 150e-6; si = 64; si1 = si;lb = 600e-6; si = 16; si1 = si;%lb = 600e-6; si = 8; si1 = si;%lb = 150e-6; si = 32; si1 = si;%lb = 150e-6; si = 128; si1 = 32;td1start = 2;%td1end = 16;thresh = .1;LoadSer = 0;SaveImagesRaw = 1;CheckThresh = 0;MakeFit = 1;PlotInterm = 0;maxfail = 5;maxiter = 20;NBS = 10;Nnodes = 100;fs = 15;cd(DataDir)%GetExpnamsfor ndir = 1:length(ExpNam)    ExpNam{ndir}    cd(ExpNam{ndir})    if exist('expno') == 0        GetExpnos    end    for nexp = 1:length(expno)        %expno(nexp)        ConvertAcqus = 'Y';    ConvertProcs = 'N';    MakeTextfile = 'N';        if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0            ConvertAcqus = 'Y';        end        res = fExpnoInfo2(ConvertAcqus,MakeTextfile,ConvertProcs,expno(nexp));        eval(['load ' num2str(expno(nexp)) '/NMRacqus'])        if any(strcmp(NMRacqus.pulprog,{'DT_moacrare2drs'})) == 1            if LoadSer                            RARE2Dtd1                                load([num2str(expno(nexp)) '/SamplingVector'])                ramp = sv.bT.ramp;                G.xa = ramp.xa*Gmax.*NMRacqus.cnst1/100;                G.xb = ramp.xb*Gmax.*NMRacqus.cnst1/100;                G.xc = ramp.xc*Gmax.*NMRacqus.cnst1/100;                G.ya = ramp.ya*Gmax.*NMRacqus.cnst2/100;                G.yb = ramp.yb*Gmax.*NMRacqus.cnst2/100;                G.yc = ramp.yc*Gmax.*NMRacqus.cnst2/100;                G.za = ramp.za*Gmax.*NMRacqus.cnst3/100;                G.zb = ramp.zb*Gmax.*NMRacqus.cnst3/100;                G.zc = ramp.zc*Gmax.*NMRacqus.cnst3/100;                %symmetry vector of b-matrix                symv.x = G.xc;                symv.y = G.yc;                symv.z = G.zc;    %             symv.x = cos(DiffRamp.bT.phi).*sin(DiffRamp.bT.theta);    %             symv.y = sin(DiffRamp.bT.phi).*sin(DiffRamp.bT.theta);    %             symv.z = cos(DiffRamp.bT.theta);                symv.norm = sqrt(symv.x.^2 + symv.y.^2 + symv.z.^2);                symv.x = symv.x./symv.norm;                symv.y = symv.y./symv.norm;                symv.z = symv.z./symv.norm;                symv.theta = acos(symv.z);                symv.phi = atan2(symv.y,symv.x);                            %gradient time-modulation                load([num2str(expno(nexp)) '/GradientModulation'])                Ndt = length(gm.c);                tau = NMRacqus.d3;                dt = tau/Ndt;                %figure(1), clf, plot((1:Ndt)',[gm.c gm.b gm.a],'-'), return                G.x = repmat(gm.a,[1, td1]).*repmat(G.xa',[Ndt, 1]) + ...                    repmat(gm.b,[1, td1]).*repmat(G.xb',[Ndt, 1]) + ...                    repmat(gm.c,[1, td1]).*repmat(G.xc',[Ndt, 1]);                G.y = repmat(gm.a,[1, td1]).*repmat(G.ya',[Ndt, 1]) + ...                    repmat(gm.b,[1, td1]).*repmat(G.yb',[Ndt, 1]) + ...                    repmat(gm.c,[1, td1]).*repmat(G.yc',[Ndt, 1]);                G.z = repmat(gm.a,[1, td1]).*repmat(G.za',[Ndt, 1]) + ...                    repmat(gm.b,[1, td1]).*repmat(G.zb',[Ndt, 1]) + ...                    repmat(gm.c,[1, td1]).*repmat(G.zc',[Ndt, 1]);                %dephasing vector F                F.x = cumsum(G.x*dt,1);                F.y = cumsum(G.y*dt,1);                F.z = cumsum(G.z*dt,1);                F.r = sqrt(F.x.^2 + F.y.^2 + F.z.^2);                %diffusion weighting matrix b                %factor 2 from the double DW blocks                bT.xx = 2*gamma^2*sum(F.x.*F.x*dt,1)';                bT.xy = 2*gamma^2*sum(F.x.*F.y*dt,1)';                bT.xz = 2*gamma^2*sum(F.x.*F.z*dt,1)';                bT.yy = 2*gamma^2*sum(F.y.*F.y*dt,1)';                bT.yz = 2*gamma^2*sum(F.y.*F.z*dt,1)';                bT.zz = 2*gamma^2*sum(F.z.*F.z*dt,1)';                b = (bT.xx + bT.yy + bT.zz);                %figure(1), clf, plot((1:td1)',[b],'o'), return                %consistency check of the b-matrix eigenvalues                G.a = sqrt(G.xa.^2 + G.ya.^2 + G.za.^2);                G.b = sqrt(G.xb.^2 + G.yb.^2 + G.zb.^2);                G.c = sqrt(G.xc.^2 + G.yc.^2 + G.zc.^2);                %figure(1), clf, plot((1:td1)',G.a,'rx',(1:td1)',G.b,'go',(1:td1)',G.c,'bs'), return                lambda.min = zeros(td1,1);                lambda.mid = zeros(td1,1);                lambda.max = zeros(td1,1);                lambda.XX = zeros(td1,1);                lambda.YY = zeros(td1,1);                lambda.ZZ = zeros(td1,1);                for ntd1 = 1:td1                    lambdas = eig([bT.xx(ntd1) bT.xy(ntd1) bT.xz(ntd1)                    bT.xy(ntd1) bT.yy(ntd1) bT.yz(ntd1)                    bT.xz(ntd1) bT.yz(ntd1) bT.zz(ntd1)]);                    [dummy,indx] = sort(lambdas,'descend');                    lambda.max(ntd1,1) = min(lambdas(indx(1)));                                            lambda.mid(ntd1,1) = min(lambdas(indx(2)));                                            lambda.min(ntd1,1) = min(lambdas(indx(3)));                                            Dlambdas = abs(lambdas-sum(lambdas)/3);                    [dummy,indx] = sort(Dlambdas,'descend');                    lambda.ZZ(ntd1,1) = lambdas(indx(1));                    lambda.XX(ntd1,1) = lambdas(indx(2));                                            lambda.YY(ntd1,1) = lambdas(indx(3));                end                lambda.trace = (lambda.XX + lambda.YY + lambda.ZZ);                lambda.iso = lambda.trace/3;                lambda.Delta = (lambda.ZZ - (lambda.YY+lambda.XX)/2)./(3*lambda.iso);                lambda.eta = (lambda.YY - lambda.XX+eps)./(2*lambda.iso.*lambda.Delta+eps);                lambda.S = lambda.min;                lambda.P = lambda.mid - lambda.min;                lambda.L = lambda.max - lambda.mid;                G.angle = acos(sqrt((2*lambda.Delta + 1)/3));                bT.N = td1;                bT.trace = lambda.trace;                bT.iso = lambda.iso;                bT.XX = lambda.XX;                bT.YY = lambda.YY;                bT.ZZ = lambda.ZZ;                bT.Delta = lambda.Delta;                bT.eta = lambda.eta;                bT.S = lambda.S;                bT.P = lambda.P;                bT.L = lambda.L;                bT.dir.x = symv.x;                bT.dir.y = symv.y;                bT.dir.z = symv.z;                bT.theta = symv.theta;                bT.phi = symv.phi;                bT.zeta = G.angle;                                                tperacq = 15*NMRacqus.d22 + sv.vdT1 + 2e-6 + 20*NMRacqus.d32 + 4e-6*NMRacqus.p11 + ...                    NMRacqus.d34 + sv.vdT2 + 2*NMRacqus.d3 + 2*NMRacqus.d4 + 1e-6*NMRacqus.p12 + ...                    NMRacqus.d35 + NMRacqus.aq + NMRacqus.d11;                indx = 1 + NMRacqus.nbl*(1:(NMRacqus.l2-1));                tsatperacq = tperacq - (.5e-6*NMRacqus.p11+3*NMRacqus.d32+2*NMRacqus.d22+NMRacqus.d11);                tperacq(indx) = tperacq(indx) + (NMRacqus.d22 + NMRacqus.d1);                t0acq = [0; cumsum(tperacq(1:(sv.N-1)))];%                 tT2W = .5e-6*NMRacqus.p11 + 7*NMRacqus.d32 + 2*NMRacqus.d22 + sv.vdT2 + NMRacqus.d34 + ...%                     2*NMRacqus.d3 + 2*NMRacqus.d4 + 1e-6*NMRacqus.p12 + NMRacqus.d35 + .5*NMRacqus.aq;                tT2W = .5e-6*NMRacqus.p11 + 7*NMRacqus.d32 + 2*NMRacqus.d22 + sv.vdT2 + NMRacqus.d34 + ...                    2*NMRacqus.d3 + 2*NMRacqus.d4 + 1e-6*NMRacqus.p12 + NMRacqus.d35 + 0*NMRacqus.aq;                tT1R = zeros(sv.N,1);                for nslice = 1:sv.Nslices                    indx = find(nslice == sv.sliceindx);                    tT1R(indx(1)) = 1/eps;                    for n = 2:sv.NRDenc                        tT1R(indx(n)) = t0acq(indx(n)) - (t0acq(indx(n-1))+tsatperacq(indx(n-1))) + 5*NMRacqus.d22 + sv.vdT1(indx(n)) + ...                            2e-6 + NMRacqus.d32 + .5e-6*NMRacqus.p11;                    end                end                       %%                          figure(2), clf                subplot(5,1,1), plot(1:sv.N,t0acq)                subplot(5,1,2), plot(1:sv.N,tperacq)                subplot(5,1,3), plot(1:sv.N,tsatperacq)                subplot(5,1,4), plot(1:sv.N,tT2W)                subplot(5,1,5), plot(1:sv.N,tT1R), set(gca,'XLim',[0 sv.N],'YLim',max(tT1R((sv.Nslices+1):end))*[0 1.1])                    %%                         if SaveImagesRaw                    ImagesRaw.S = Itd1;                    ImagesRaw.r = r;                    ImagesRaw.proc.lb = lb;                    ImagesRaw.bT = bT;                    ImagesRaw.tT2W = tT2W;                    ImagesRaw.tT1R = tT1R;                    eval(['save ' num2str(expno(nexp)) '/ImagesRaw ImagesRaw'])                end            else                eval(['load ' num2str(expno(nexp)) '/ImagesRaw'])                load([num2str(expno(nexp)) '/SamplingVector'])                load([num2str(expno(nexp)) '/GradientModulation'])                Itd1 = ImagesRaw.S;                r = ImagesRaw.r;                lb = ImagesRaw.proc.lb;                bT = ImagesRaw.bT;                tT2W = ImagesRaw.tT2W;                tT1R = ImagesRaw.tT1R;                nudim.i = numel(r.i);                nudim.j = numel(r.j);                td1 = sv.N;            end%%            Imax = abs(Itd1);            Imax(nudim.i/2+(0:2),nudim.j/2+(0:2),2:td1) = 0;            Imax(:,:,1:(td1start-1)) = 0;            Imax = max(reshape(Imax,numel(Imax),1));                            rmax.i = max(r.i)*1e3;            rmax.j = max(r.j)*1e3;%%                         %figure(1), clf, plot((1:td1)',bT.XX,'r-',(1:td1)',bT.YY,'g--',(1:td1)',bT.ZZ,'b-',(1:td1)',bT.trace,'k--'), return            %figure(1), clf, plot((1:td1)',bT.Delta,'r-',(1:td1)',bT.eta,'g-'), return%             figure(3), clf%             subplot(4,2,1)%             plot(1:td1,bT.trace)%             subplot(4,2,3)%             plot(1:td1,bT.Delta)%             subplot(2,2,2)%             plot(1:length(gm.a),gm.a,'r-',1:length(gm.a),gm.b,'g-',...%                 1:length(gm.a),gm.c,'k-')%             subplot(2,2,3)%             plot(1:td1,bT.dir.x,'r.',1:td1,bT.dir.y,'g.',1:td1,bT.dir.z,'b.')% %             [X,Y] = fSchmidt(bT.dir.x,bT.dir.y,bT.dir.z);%             latitude.theta = pi/180*[30:30:150 179];%             latitude.phi = linspace(0,2*pi,100);%             [latitude.phi,latitude.theta] = ndgrid(latitude.phi,latitude.theta);%             latitude.z = cos(latitude.theta);%             latitude.x = sin(latitude.theta).*cos(latitude.phi);%             latitude.y = sin(latitude.theta).*sin(latitude.phi);%             [latitude.X,latitude.Y] = fSchmidt(latitude.x,latitude.y,latitude.z);%             longitude.theta = pi/180*linspace(30,180,100);%             longitude.phi = pi/180*[30:30:360];%             [longitude.theta,longitude.phi] = ndgrid(longitude.theta,longitude.phi);%             longitude.z = cos(longitude.theta);%             longitude.x = sin(longitude.theta).*cos(longitude.phi);%             longitude.y = sin(longitude.theta).*sin(longitude.phi);%             [longitude.X,longitude.Y] = fSchmidt(longitude.x,longitude.y,longitude.z);% %             subplot(2,2,4)        %             plot(X,Y,'k.')%             hold on%             plot(latitude.X,latitude.Y,'b-')%             plot(longitude.X,longitude.Y,'b-')%             axis equal %%                                   %             figure(4), clf%             colormap(hot(256));%             clim = [0 1];%             %             td1dim.x = ceil(11/8.5*sqrt(td1));%             td1dim.y = ceil(8.5/11*sqrt(td1));%             if any(strcmp(NMRacqus.pulprog,{'DT_dtirare2d','DT_dpgsedtirare2d'})) == 1                        %                 td1dim.x = 8;%                 td1dim.y = ceil(td1/td1dim.x);%             end%             if any(strcmp(NMRacqus.pulprog,{'DT_rare2dms','DT_moacrare2dms'})) == 1                        %                 td1dim.x = NMRacqus.nbl;%                 td1dim.y = ceil(td1/td1dim.x);%             end%             td1dim.x = 15;%             td1dim.y = 10;%             width = 1/td1dim.x;%             height = 1/td1dim.y;%             ntd1 = 1;%             for ycount = 1:td1dim.y;%                 for xcount = 1:td1dim.x;%                     if ntd1 <= td1;%                         axes('position',[(xcount-1)*width 1-ycount*height width height])%                         Iplot = squeeze(abs(Itd1(:,:,ntd1)))/Imax;%                         Iplot = Iplot/max(max(Iplot));% %                         if tdim.i>tdim.j% %                             imagesc(r.j*1e3,r.i*1e3,Iplot,clim)% %                         else%                             imagesc(r.i*1e3,r.j*1e3,Iplot',clim)% %                         end%                         set(gca,'YDir','normal')%                         axis equal, axis tight, axis off%                     end%                     ntd1 = ntd1+1;%                 end%             end            %return            %             figure(1), clf%             Iplot = abs(squeeze(sum(Itd1(:,:,1:sv.Nslices),3)));%             colormap('gray')%             Iplot = Iplot/max(max(Iplot));%             Iplot(Iplot>thresh) = 1;%             clim = [0 1];%             imagesc(r.i*1e3,r.j*1e3,Iplot',clim)%             set(gca,'YDir','normal')%             axis equal, axis tight% %             if CheckThresh%                 return%             end%%            if MakeFit       %%                Dmin = 1e-12; Dmax = 3e-9; R1min = .1; R1max = 2; R2min = 2; R2max = 50;                                rng(1,'twister')                Npoints = sv.NRDenc;                NpointsBS = Npoints;                %NpointsBS = min([512; Npoints]);                nudim_read = nudim.i;                nudim_phase = nudim.j;                nudim_slice = sv.Nslices;                                nslice = 8;                sliceindx = find(nslice == sv.sliceindx);                Sslice = ImagesRaw.S(:,:,sliceindx);                DRenc.tT1R = tT1R(sliceindx);                DRenc.tT1W = zeros(sv.NRDenc,1);                DRenc.tT2W = tT2W(sliceindx);                DRenc.bT.trace = bT.trace(sliceindx);                DRenc.bT.Delta = bT.Delta(sliceindx);                DRenc.bT.theta = bT.theta(sliceindx);                DRenc.bT.phi = bT.phi(sliceindx);                td1dim.x = ceil(11/8*sqrt(sv.NRDenc));                td1dim.y = ceil(8/11*sqrt(sv.NRDenc));                td1dim.x = 10;                td1dim.y = 8;                width = 1/td1dim.x;                height = 1/td1dim.y;                ntd1 = 1;                figure(1), clf                nread = 13; nphase = 6;                nread = 6; nphase = 10;                %nread = 10; nphase = 8;                for ycount = 1:td1dim.y;                    for xcount = 1:td1dim.x;                        if ntd1 <= sv.NRDenc;                            axes('position',[(xcount-1)*width 1-ycount*height width height])                            Iplot = squeeze(abs(Sslice(:,:,ntd1)))/Imax;                            %Iplot = Iplot/max(max(Iplot));    %                         if tdim.i>tdim.j    %                             imagesc(r.j*1e3,r.i*1e3,Iplot,clim)    %                         else                                imagesc(r.i*1e3,r.j*1e3,Iplot',[0 1])                                hold on, plot(r.i(nread)*1e3,r.j(nphase)*1e3,'rx')    %                         end                            set(gca,'YDir','normal')                            axis equal, axis tight, axis off                        end                        ntd1 = ntd1+1;                    end                end                %%                Signal = abs(squeeze(Sslice(nread,nphase,:)));                figure(2), clf                subplot(7,1,1), plot(1:sv.NRDenc,DRenc.tT1R), set(gca,'YLim',max(tT1R((sv.Nslices+1):end))*[0 1.1])                             subplot(7,1,2), plot(1:sv.NRDenc,DRenc.tT2W)                subplot(7,1,3), plot(1:sv.NRDenc,DRenc.bT.trace)                subplot(7,1,4), plot(1:sv.NRDenc,DRenc.bT.Delta)                subplot(7,1,5), plot(1:sv.NRDenc,DRenc.bT.theta)                subplot(7,1,6), plot(1:sv.NRDenc,DRenc.bT.phi)                subplot(7,1,7), plot(1:sv.NRDenc,Signal)return%%                                ImagesMOAC_w = zeros(nudim_read,nudim_phase,nudim_slice,Nnodes,NBS);                ImagesMOAC_iso = zeros(nudim_read,nudim_phase,nudim_slice,Nnodes,NBS);                ImagesMOAC_Delta = zeros(nudim_read,nudim_phase,nudim_slice,Nnodes,NBS);                ImagesMOAC_theta = zeros(nudim_read,nudim_phase,nudim_slice,Nnodes,NBS);                ImagesMOAC_phi = zeros(nudim_read,nudim_phase,nudim_slice,Nnodes,NBS);                ImagesMOAC_R1 = zeros(nudim_read,nudim_phase,nudim_slice,Nnodes,NBS);                ImagesMOAC_R2 = zeros(nudim_read,nudim_phase,nudim_slice,Nnodes,NBS);                ImagesMOAC_chisq = zeros(nudim_read,nudim_phase,nudim_slice,NBS);                ImagesMOAC_nnonzero = zeros(nudim_read,nudim_phase,nudim_slice,NBS);                options = optimset('MaxFunEvals',1e4,'Display','off');                p =  TimedProgressBar( nudim_slice*nudim_slice*nudim_read, 10, ...                'Computing. Remaining time: ', ', Completed: ', 'Concluded in ' );                %for nslice = 1:nudim_slice;                for nslice = 8;                    sliceindx = find(nslice == sv.sliceindx);                    Sslice = ImagesRaw.S(:,:,sliceindx);                    DRenc.tT1R = tT1R(sliceindx);                    DRenc.tT1W = zeros(sv.NRDenc,1);                    DRenc.tT2W = tT2W(sliceindx);                    DRenc.bT.trace = bT.trace(sliceindx);                    DRenc.bT.Delta = bT.Delta(sliceindx);                    DRenc.bT.theta = bT.theta(sliceindx);                    DRenc.bT.phi = bT.phi(sliceindx);                    for nphase = 1:nudim_phase                            for nread = 1:nudim_read    %%                            %nread = 13; nphase = 6;                            %nread = 6; nphase = 10;                            %nread = 10; nphase = 8;                            %[nread nphase]                            Signal = abs(squeeze(Sslice(nread,nphase,:)));                            S_vector = Signal;                            S0 = max(S_vector);    %                         Diso = 5e-10; DDelta = .85;    %                         bDDelDel = bT.trace.*Diso.*bT.Delta.*DDelta;    %                         S_powder = S0*exp(-bT.trace.*Diso).*exp(bDDelDel).*...    %                         sqrt(pi)/2.*real(gammainc(3*bDDelDel,1/2)./sqrt(3*bDDelDel)).*...    %                         exp(-bT.TR1w.*R1).*(1 - exp(-bT.TR.*R1)).*exp(-bT.TE.*R2);    %                         P_powder = .95;    %                         S_vector = (1-P_powder)*S_vector + P_powder * S_powder;                            %figure(1), clf, semilogy(DRenc.bT.trace,S_vector,'o'), return                            %figure(1), clf, semilogy(DRenc.tT2W,S_vector,'o',DRenc.tT2W,7*max(S_vector)*exp(-45*DRenc.tT2W),'x'), return                            %figure(1), clf, semilogy(DRenc.tT2W,S_vector,'o',DRenc.tT2W,1*max(S_vector)*exp(-2*DRenc.tT2W),'x'), return                            chisq_v = zeros(NBS,1);                            niter_v = zeros(NBS,1);                            nfail_v = zeros(NBS,1);                            nnonzero_v = zeros(NBS,1);                            uDTrec_iso_a = zeros(Nnodes,NBS);                            uDTrec_Delta_a = zeros(Nnodes,NBS);                            uDTrec_theta_a = zeros(Nnodes,NBS);                            uDTrec_phi_a = zeros(Nnodes,NBS);                            uDTrec_R1_a = zeros(Nnodes,NBS);                            uDTrec_R2_a = zeros(Nnodes,NBS);                            uDTrec_w_a = zeros(Nnodes,NBS);                            for nBS = 1:NBS                            %parfor nBS = 1:NBS                            %%                                    %indx_BS = (1:Npoints)';                                indx_BS = td1start - 1 + sort(ceil((Npoints-td1start+1)*rand(NpointsBS,1)));                                %figure(1), clf, plot(1:NpointsBS,indx_BS,'-'), return                                indx_BS_array = repmat(indx_BS',[NpointsBS 1]);                                bT_trace = repmat(DRenc.bT.trace,1,Nnodes);                                bT_Delta = repmat(DRenc.bT.Delta,1,Nnodes);                                bT_theta = repmat(DRenc.bT.theta,1,Nnodes);                                bT_phi = repmat(DRenc.bT.phi,1,Nnodes);                                t_T1R = repmat(DRenc.tT1R,1,Nnodes);                                t_T1W = repmat(DRenc.tT1W,1,Nnodes);                                t_T2W = repmat(DRenc.tT2W,1,Nnodes);                                %Nnodes = Nnodes_v(nBS);                                chisq0 = S0^2;                                chisq = .99*S0^2;                                niter2 = 0; niter = 0;                                uDTrec_iso_v = zeros(Nnodes,1);                                uDTrec_Delta_v = zeros(Nnodes,1);                                uDTrec_theta_v = zeros(Nnodes,1);                                uDTrec_phi_v = zeros(Nnodes,1);                                uDTrec_R1_v = zeros(Nnodes,1);                                uDTrec_R2_v = zeros(Nnodes,1);                                uDTrec_w_v = zeros(Nnodes,1);                    %            while all([sqrt(chisq_temp) > chilim niter2 < maxiter2])                                niter = 0;                                uDTrec_w_v = 0*uDTrec_w_v;                                nfail = 0;                                niter = 0;                                while all([nfail < maxfail niter < maxiter])                                    uDTrec_par_vtemp = Dmin*(Dmax/Dmin).^rand(Nnodes,1);                                    uDTrec_perp_vtemp = Dmin*(Dmax/Dmin).^rand(Nnodes,1);                                    %indx_iso = round(Nnodes*.9):Nnodes;                                    %indx_iso = ceil(rand(ceil(.1*Nnodes),1)*Nnodes);                                    %uDTrec_perp_vtemp(indx_iso) = uDTrec_par_vtemp(indx_iso);                                    uDTrec_theta_vtemp = acos(2*rand(Nnodes,1)-1);                                    uDTrec_phi_vtemp = 2*pi*rand(Nnodes,1);                                    %indx_powder = ceil(rand(ceil(.3*Nnodes),1)*Nnodes);                                    %uDTrec_theta_vtemp(indx_powder) = NaN;                                    uDTrec_R1_vtemp = R1min*(R1max/R1min).^rand(Nnodes,1);                                    uDTrec_R2_vtemp = R2min*(R2max/R2min).^rand(Nnodes,1);%                                     uDTrec_par_vtemp =  Dpar*ones(Nnodes,1);%                                     uDTrec_perp_vtemp =  Dperp*ones(Nnodes,1);%                                     uDTrec_R1_vtemp = R1*ones(Nnodes,1);%                                     uDTrec_R2_vtemp = R2*ones(Nnodes,1);%                                     uDTrec_par_vtemp = max([uDTrec_par_vtemp uDTrec_perp_vtemp],[],2);                                    [uDTrec_w_v,indx] = sort(uDTrec_w_v,'descend');                                    uDTrec_iso_v = uDTrec_iso_v(indx);                                    uDTrec_Delta_v = uDTrec_Delta_v(indx);                                    uDTrec_theta_v = uDTrec_theta_v(indx);                                    uDTrec_phi_v = uDTrec_phi_v(indx);                                    uDTrec_R1_v = uDTrec_R1_v(indx);                                    uDTrec_R2_v = uDTrec_R2_v(indx);                                    uDTrec_par_v = uDTrec_iso_v.*(1 + 2*uDTrec_Delta_v);                                    uDTrec_perp_v = uDTrec_iso_v.*(1 - uDTrec_Delta_v);                                    nnonzero = sum(uDTrec_w_v>0);                                    indx = 1:nnonzero;                                    uDTrec_par_vtemp(indx) = uDTrec_par_v(indx);                                    uDTrec_perp_vtemp(indx) = uDTrec_perp_v(indx);                                    uDTrec_theta_vtemp(indx) = uDTrec_theta_v(indx);                                    uDTrec_phi_vtemp(indx) = uDTrec_phi_v(indx);                                    uDTrec_R1_vtemp(indx) = uDTrec_R1_v(indx);                                    uDTrec_R2_vtemp(indx) = uDTrec_R2_v(indx);                                    %indx = 1:min([nnonzero floor(Nnodes/10)]);                                    %indx2 = nnonzero + (1:numel(indx));                                    indx2 = indx+min([nnonzero (Nnodes-max(indx))]);                                    uDTrec_par_vtemp(indx2) = uDTrec_par_v(indx).*(1+.1*randn(nnonzero,1));                                    uDTrec_perp_vtemp(indx2) = uDTrec_perp_v(indx).*(1+.1*randn(nnonzero,1));                                    uDTrec_theta_vtemp(indx2) = uDTrec_theta_v(indx) + .01*2*pi*randn(nnonzero,1);                                    uDTrec_phi_vtemp(indx2) = uDTrec_phi_v(indx) + .01*2*pi*randn(nnonzero,1);                                    uDTrec_R1_vtemp(indx2) = uDTrec_R1_v(indx).*(1+.1*randn(nnonzero,1));                                    uDTrec_R2_vtemp(indx2) = uDTrec_R2_v(indx).*(1+.1*randn(nnonzero,1));                                    uDTrec_par_vtemp(uDTrec_par_vtemp>Dmax) = Dmax;                                    uDTrec_par_vtemp(uDTrec_par_vtemp<Dmin) = Dmin;                                    uDTrec_perp_vtemp(uDTrec_perp_vtemp>Dmax) = Dmax;                                    uDTrec_perp_vtemp(uDTrec_perp_vtemp<Dmin) = Dmin;                                    uDTrec_R1_vtemp(uDTrec_R1_vtemp>R1max) = R1max;                                    uDTrec_R1_vtemp(uDTrec_R1_vtemp<R1min) = R1min;                                    uDTrec_R2_vtemp(uDTrec_R2_vtemp>R2max) = R2max;                                    uDTrec_R2_vtemp(uDTrec_R2_vtemp<R2min) = R2min;                                                                        uDTrec_iso_vtemp = (uDTrec_par_vtemp + 2*uDTrec_perp_vtemp)/3;                                    uDTrec_Delta_vtemp = (uDTrec_par_vtemp - uDTrec_perp_vtemp)/3./uDTrec_iso_vtemp;                                    uDTrec_iso = repmat(uDTrec_iso_vtemp',Npoints,1);                                    uDTrec_Delta = repmat(uDTrec_Delta_vtemp',Npoints,1);                                    uDTrec_theta = repmat(uDTrec_theta_vtemp',Npoints,1);                                    uDTrec_phi = repmat(uDTrec_phi_vtemp',Npoints,1);                                    uDTrec_R1 = repmat(uDTrec_R1_vtemp',Npoints,1);                                    uDTrec_R2 = repmat(uDTrec_R2_vtemp',Npoints,1);                                    cosbeta = cos(bT_theta).*cos(uDTrec_theta) + sin(bT_theta).*sin(uDTrec_theta).*cos(bT_phi - uDTrec_phi);                                    P2cosbeta = (3*cosbeta.^2 - 1)/2;                                    Deff = uDTrec_iso.*(1 + 2*bT_Delta.*uDTrec_Delta.*P2cosbeta);                                    Kernel = exp(-bT_trace.*Deff).*exp(-t_T1W.*uDTrec_R1)...                                        .*(1 - exp(-t_T1R.*uDTrec_R1)).*exp(-t_T2W.*uDTrec_R2);    %                                 indx_powder = find(isnan(uDTrec_theta_vtemp));    %                                 bDDelDel = bT_trace.*uDTrec_iso.*bT_Delta.*uDTrec_Delta;    %     %                                 Kernel_powder = exp(-bT_trace.*uDTrec_iso).*exp(bDDelDel).*...    %                                     sqrt(pi)/2.*real(gammainc(3*bDDelDel,1/2)./sqrt(3*bDDelDel)).*...    %                                     exp(-tT1W.*uDTrec_R1).*(1 - exp(-tT1R.*uDTrec_R1)).*exp(-tT2W.*uDTrec_R2);    %     %                                 indx = bDDelDel == 0;    %                                 Kernel_powder(indx) = exp(-bT_trace(indx).*uDTrec_iso(indx));    %                                 Kernel_powder(bT_trace == 0) = 1;    %                                 Kernel_powder(bDDelDel < -10) = 0;    %                                 indx = isnan(Kernel_powder);    %                                 Kernel_powder(indx) = 0;    %                                 indx = isinf(Kernel_powder);    %                                 Kernel_powder(indx) = 0;    %                                     %                                 Kernel(:,indx_powder) = Kernel_powder(:,indx_powder);                                     %figure(1), clf, plot(t_T2W,Kernel,'-'), return                                    S_vector_temp = S_vector(indx_BS,1);                                    Kernel_temp = Kernel(indx_BS,:);                    %                     [S_vector_temp,indx] = sort(S_vector_temp,1,'descend');                    %                     Kernel_temp = Kernel_temp(indx,:);                                    uDTrec_w_vtemp = S0*lsqnonneg(Kernel_temp,S_vector_temp/S0);                                    Scalc = Kernel_temp*uDTrec_w_vtemp;                                    chisq = sum((S_vector_temp-Scalc).^2,1)/Npoints;    %     %                                 figure(2), clf    %                                 plot(1:NpointsBS,Scalc,'k-',1:NpointsBS,S_vector_temp,'bo',[0; NpointsBS],sqrt(chisq)*[1; 1],'k--')    %                                 title([num2str(niter) ' ' num2str(nfail) ' ' num2str(nnonzero) ' ' num2str(chisq,3) ' ' num2str(chisq0,3)])    %                                 pause(.1)                                    niter = niter+1;                                    if chisq > chisq0*(1 - 1e-3)                                        nfail = nfail+1;                                    end                                    if chisq < chisq0*(1 - 0e-2)                                        uDTrec_iso_v = uDTrec_iso_vtemp;                                        uDTrec_Delta_v = uDTrec_Delta_vtemp;                                        uDTrec_theta_v = uDTrec_theta_vtemp;                                        uDTrec_phi_v = uDTrec_phi_vtemp;                                        uDTrec_R1_v = uDTrec_R1_vtemp;                                        uDTrec_R2_v = uDTrec_R2_vtemp;                                        uDTrec_w_v = uDTrec_w_vtemp;                                        chisq0 = chisq;                                    end                                end                                [uDTrec_w_v,indx] = sort(uDTrec_w_v,'descend');                                uDTrec_iso_v = uDTrec_iso_v(indx);                                uDTrec_Delta_v = uDTrec_Delta_v(indx);                                uDTrec_theta_v = uDTrec_theta_v(indx);                                uDTrec_phi_v = uDTrec_phi_v(indx);                                 uDTrec_R1_v = uDTrec_R1_v(indx);                                 uDTrec_R2_v = uDTrec_R2_v(indx);                                 nnonzero = sum(uDTrec_w_v>0);                                uDTrec_iso_a(:,nBS) = uDTrec_iso_v;                                uDTrec_Delta_a(:,nBS) = uDTrec_Delta_v;                                uDTrec_theta_a(:,nBS) = uDTrec_theta_v;                                uDTrec_phi_a(:,nBS) = uDTrec_phi_v;                                uDTrec_R1_a(:,nBS) = uDTrec_R1_v;                                uDTrec_R2_a(:,nBS) = uDTrec_R2_v;                                uDTrec_w_a(:,nBS) = uDTrec_w_v;                                chisq_v(nBS,1) = chisq;                                niter_v(nBS,1) = niter;                                nfail_v(nBS,1) = nfail;                                nnonzero_v(nBS,1) = nnonzero;                    %                 figure(2), clf                    %                 plot(1:Npoints,S_vector_temp,'o',1:Npoints,Scalc,'-')                    %                 pause(.1)                            end                    %             figure(1), clf                    %             subplot(2,2,1)                    %             plot(1:NBS,sqrt(chisq_v),'o')                    %             ylabel('rms chi'), xlabel('nBS')                    %             subplot(2,2,2)                    %             hist(sqrt(chisq_v))                    %             xlabel('chisq')                    %             subplot(2,2,3)                    %             hist(niter_v)                    %             xlabel('niter')                    %             subplot(2,2,4)                    %             hist(nnonzero_v)                    %             xlabel('nnonzero')                            ImagesMOAC_w(nread,nphase,nslice,:,:) = uDTrec_w_a;                            ImagesMOAC_iso(nread,nphase,nslice,:,:) = uDTrec_iso_a;                            ImagesMOAC_Delta(nread,nphase,nslice,:,:) = uDTrec_Delta_a;                            ImagesMOAC_theta(nread,nphase,nslice,:,:) = uDTrec_theta_a;                            ImagesMOAC_phi(nread,nphase,nslice,:,:) = uDTrec_phi_a;                            ImagesMOAC_R1(nread,nphase,nslice,:,:) = uDTrec_R1_a;                            ImagesMOAC_R2(nread,nphase,nslice,:,:) = uDTrec_R2_a;                            ImagesMOAC_chisq(nread,nphase,nslice,:) = chisq_v;                            ImagesMOAC_nnonzero(nread,nphase,nslice,:) = nnonzero_v;                            p.progress; %Counter for progress report                        end                    end                end                p.stop;          %%                ImagesMOAC_a.w = ImagesMOAC_w;                ImagesMOAC_a.iso = ImagesMOAC_iso;                ImagesMOAC_a.Delta = ImagesMOAC_Delta;                ImagesMOAC_a.theta = ImagesMOAC_theta;                ImagesMOAC_a.phi = ImagesMOAC_phi;                ImagesMOAC_a.R1 = ImagesMOAC_R1;                ImagesMOAC_a.R2 = ImagesMOAC_R2;                ImagesMOAC_a.chisq = ImagesMOAC_chisq;                ImagesMOAC_a.nnonzero = ImagesMOAC_nnonzero;                ReconDat = struct('NBS',NBS,'Nnodes',Nnodes,'Dmin',Dmin,'Dmax',Dmax,...                    'R1min',R1min,'R1max',R1max,'R2min',R2min,'R2max',R2max);                save([num2str(expno(nexp)) '/ImagesMOAC_a'],'ImagesMOAC_a','ReconDat','-v7.3')                   end        end    end        cd ..endcd(wd)