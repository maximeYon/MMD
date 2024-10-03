clear all

wd = cd;

DataDir = '/Users/daniel/NMRdata/AVII500/DT/';
%DataDir = '/Users/daniel/NMRdata/AVII200/';
%DataDir = '/opt/topspin2/data/DT/nmr';
%DataDir = '/opt/topspin/data/DT/nmr';
%DataDir = '/Users/daniel/Dropbox';
%DataDir = '/Users/daniel/Documents/Spaces/Presentations';

% ExpNam = {'RandSampDR1R2'}; expno = [30 31 33 49:55 66 68:71 93 96 99 102 105 107 109 111:2:119 128 130];
% expno = [141];
ExpNam = {'AOToct_temp10'}; expno = [17:147]; expno = [157:10:217]; expno = 37;
%ExpNam = {'AOToct8'}; expno = 16;
%ExpNam = {'AOToct9'}; expno = 20;
%ExpNam = {'C14E5_6'};
%ExpNam = {'C14E5_7'}; %expno = [16 19]; expno = 22;
%ExpNam = {'AOToct_Eq6'}; expno = 20:10:240;
%ExpNam = {'AOToct_Eq7'}; expno = 20:10:240;

AutoPhase = 1;
FindPeaks = 1;
CheckBasline = 0;
CheckPeaks = 0;
PlotInterm = 0;
MakeMCfit = 1;

%optional process parameters
thresh = .1;
td1start = 1;
lb = 10; %si = 4*1024;
% basl = [.05 .15 .4 .6 .85 .95]; %skip = [0 .61 .67 1];
% AutoPhase = 1; %phc0 = 30; phc1 = 0; pivot = .5;

Imin = .8e-2; %min intensity for I vs b fig
signal = 'area'; %signal = 'intensity';

Dmin = 1e-11; Dmax = 5e-9; R1min = .1; R1max = 10; R2min = .5; R2max = 500;
%Dmin = 5e-11; Dmax = 5e-9; R1min = .1; R1max = 1; R2min = 5; R2max = 500;

maxfail = 50;
maxiter = 500;
NBS = 1000;
Nnodes = 500;

cd(DataDir)

cd(ExpNam{1})
if exist('expno') == 0
    GetExpnos
end
for nexp = 1:length(expno)
    ConvertAcqus = 'N';    ConvertProcs = 'N';    MakeTextfile = 'N';
    if exist([num2str(expno(nexp)) '/NMRacqus.mat']) == 0
        ConvertAcqus = 'Y';
    end
    res = fExpnoInfo2(ConvertAcqus,MakeTextfile,ConvertProcs,expno(nexp));
    eval(['load ' num2str(expno(nexp)) '/NMRacqus'])

    if any(strcmp(NMRacqus.pulprog,{'DT_tbppgsteT1T2'})) == 1
        
        %FFT, phase correction, baseline correction
        Spec2Dtd1
        
        %Peak picking
        PeakPick        
        save([num2str(expno(nexp)) '/PP'],'PP')
        %Signal = sum(Apeak,2);
        Signal = Apeak(:,Npeaks);
        
        %Calculate b-values
        
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
        
        %timing variables
        epsilon = NMRacqus.d2;
        tau = 2*NMRacqus.d4 + 1e-6*NMRacqus.p2;
        delta = 2*(NMRacqus.d2 + NMRacqus.d3);
        Delta = 2*NMRacqus.d2 + 2*NMRacqus.d3 + 2*NMRacqus.d4 + 2*NMRacqus.d6 + NMRacqus.d5;            
        tT2 = 12*(NMRacqus.d6 + 2*NMRacqus.d2 + NMRacqus.d3 + NMRacqus.d4);
        tT1 = 3*NMRacqus.d5;
        if NMRacqus.l11
            Delta = Delta + NMRacqus.d46 + 2*NMRacqus.d42 + NMRacqus.d43;
            tT1 = tT1 + 3*(NMRacqus.d46 + 2*NMRacqus.d42 + NMRacqus.d43);
        end
        tdiff = Delta - delta/3 - tau/2 - epsilon/2 - epsilon^2/6/delta + epsilon^3/15/delta^2;
        
        bmax = 3*gamma^2*(Gmax.*NMRacqus.cnst1/100)^2*delta^2*tdiff;

        eval(['load ' num2str(expno(nexp)) '/DiffRamp'])
        
        bT = DiffRamp.bT;
        bT.trace = bT.trace*bmax;
        bT.TE = bT.TE + 0*tT2;
        bT.TR1w = 0*tT1.*ones(bT.N,1);
        
%         [Signal,indx] = sort(Signal,1,'descend');
%         
%         bT.trace = bT.trace(indx);
%         bT.Delta = bT.Delta(indx);
%         bT.theta = bT.theta(indx);
%         bT.phi = bT.phi(indx);
%         bT.TE = bT.TE(indx);
%         bT.TR = bT.TR(indx);
%         bT.TR1w = bT.TR1w(indx);

        figure(1), clf
        subplot(2,2,1), semilogy(1:bT.N,Signal,'-'), axis tight
        subplot(2,2,2), plot(1:bT.N,bT.trace,'-'), axis tight
        subplot(2,2,3), plot(1:bT.N,bT.TE,'-'), axis tight
        subplot(2,2,4), plot(1:bT.N,bT.TR,'-'), axis tight
        pause(.1)
        %return

        if MakeMCfit

            rng(1,'twister')
            Npoints = bT.N;
            NpointsBS = Npoints;
            %NpointsBS = min([512; Npoints]);

            S_vector = Signal;
            S0 = max(S_vector);
            %figure(1), clf, plot(bT.trace,S_vector,'-'), return

            chisq_v = zeros(NBS,1);
            niter_v = zeros(NBS,1);
            nfail_v = zeros(NBS,1);
            nnonzero_v = zeros(NBS,1);

            uDTrec_iso_a = zeros(Nnodes,NBS);
            uDTrec_Delta_a = zeros(Nnodes,NBS);
            uDTrec_theta_a = zeros(Nnodes,NBS);
            uDTrec_phi_a = zeros(Nnodes,NBS);
            uDTrec_R1_a = zeros(Nnodes,NBS);
            uDTrec_R2_a = zeros(Nnodes,NBS);
            uDTrec_w_a = zeros(Nnodes,NBS);

            %Parallel fitting
            p =  TimedProgressBar( NBS, 10, ...
            'Computing. Remaining time: ', ', Completed: ', 'Concluded in ' );
            parfor nBS = 1:NBS
            %for nBS = 1:NBS
            %%    
                %indx_BS = (1:Npoints)';
                indx_BS = td1start - 1 + sort(ceil((Npoints-td1start+1)*rand(NpointsBS,1)));
                %figure(1), clf, plot(1:NpointsBS,indx_BS,'-'), return
                indx_BS_array = repmat(indx_BS',[NpointsBS 1]);

                bT_trace = repmat(bT.trace,1,Nnodes);
                bT_Delta = repmat(bT.Delta,1,Nnodes);
                bT_theta = repmat(bT.theta,1,Nnodes);
                bT_phi = repmat(bT.phi,1,Nnodes);
                bT_TR = repmat(bT.TR,1,Nnodes);
                bT_TR1w = repmat(bT.TR1w,1,Nnodes);
                bT_TE = repmat(bT.TE,1,Nnodes);

                %Nnodes = Nnodes_v(nBS);
                chisq0 = S0^2;
                chisq = .99*S0^2;
                niter2 = 0; niter = 0;

                uDTrec_iso_v = zeros(Nnodes,1);
                uDTrec_Delta_v = zeros(Nnodes,1);
                uDTrec_theta_v = zeros(Nnodes,1);
                uDTrec_phi_v = zeros(Nnodes,1);
                uDTrec_R1_v = zeros(Nnodes,1);
                uDTrec_R2_v = zeros(Nnodes,1);
                uDTrec_w_v = zeros(Nnodes,1);

    %            while all([sqrt(chisq_temp) > chilim niter2 < maxiter2])

                    niter = 0;
                    uDTrec_w_v = 0*uDTrec_w_v;

                    nfail = 0;
                    niter = 0;
                    while all([nfail < maxfail niter < maxiter])
                        uDTrec_par_vtemp = Dmin*(Dmax/Dmin).^rand(Nnodes,1);
                        uDTrec_perp_vtemp = Dmin*(Dmax/Dmin).^rand(Nnodes,1);
                        indx_iso = round(Nnodes*.9):Nnodes;
                        uDTrec_perp_vtemp(indx_iso) = uDTrec_par_vtemp(indx_iso);
                        uDTrec_theta_vtemp = acos(2*rand(Nnodes,1)-1);
                        uDTrec_phi_vtemp = 2*pi*rand(Nnodes,1);
                        uDTrec_R1_vtemp = R1min*(R1max/R1min).^rand(Nnodes,1);
                        uDTrec_R2_vtemp = R2min*(R2max/R2min).^rand(Nnodes,1);

                        [uDTrec_w_v,indx] = sort(uDTrec_w_v,'descend');
                        uDTrec_iso_v = uDTrec_iso_v(indx);
                        uDTrec_Delta_v = uDTrec_Delta_v(indx);
                        uDTrec_theta_v = uDTrec_theta_v(indx);
                        uDTrec_phi_v = uDTrec_phi_v(indx);
                        uDTrec_R1_v = uDTrec_R1_v(indx);
                        uDTrec_R2_v = uDTrec_R2_v(indx);
                        uDTrec_par_v = uDTrec_iso_v.*(1 + 2*uDTrec_Delta_v);
                        uDTrec_perp_v = uDTrec_iso_v.*(1 - uDTrec_Delta_v);

                        nnonzero = sum(uDTrec_w_v>0);
                        indx = 1:nnonzero;
                        uDTrec_par_vtemp(indx) = uDTrec_par_v(indx);
                        uDTrec_perp_vtemp(indx) = uDTrec_perp_v(indx);
                        uDTrec_theta_vtemp(indx) = uDTrec_theta_v(indx);
                        uDTrec_phi_vtemp(indx) = uDTrec_phi_v(indx);
                        uDTrec_R1_vtemp(indx) = uDTrec_R1_v(indx);
                        uDTrec_R2_vtemp(indx) = uDTrec_R2_v(indx);

                        indx2 = indx+nnonzero;
                        uDTrec_par_vtemp(indx2) = uDTrec_par_v(indx).*(1+.1*randn(nnonzero,1));
                        uDTrec_perp_vtemp(indx2) = uDTrec_perp_v(indx).*(1+.1*randn(nnonzero,1));
                        uDTrec_theta_vtemp(indx2) = uDTrec_theta_v(indx) + .01*2*pi*randn(nnonzero,1);
                        uDTrec_phi_vtemp(indx2) = uDTrec_phi_v(indx) + .01*2*pi*randn(nnonzero,1);
                        uDTrec_R1_vtemp(indx2) = uDTrec_R1_v(indx).*(1+.1*randn(nnonzero,1));
                        uDTrec_R2_vtemp(indx2) = uDTrec_R2_v(indx).*(1+.1*randn(nnonzero,1));

                        uDTrec_iso_vtemp = (uDTrec_par_vtemp + 2*uDTrec_perp_vtemp)/3;
                        uDTrec_Delta_vtemp = (uDTrec_par_vtemp - uDTrec_perp_vtemp)/3./uDTrec_iso_vtemp;

                        uDTrec_iso = repmat(uDTrec_iso_vtemp',Npoints,1);
                        uDTrec_Delta = repmat(uDTrec_Delta_vtemp',Npoints,1);
                        uDTrec_theta = repmat(uDTrec_theta_vtemp',Npoints,1);
                        uDTrec_phi = repmat(uDTrec_phi_vtemp',Npoints,1);
                        uDTrec_R1 = repmat(uDTrec_R1_vtemp',Npoints,1);
                        uDTrec_R2 = repmat(uDTrec_R2_vtemp',Npoints,1);

                        cosbeta = cos(bT_theta).*cos(uDTrec_theta) + sin(bT_theta).*sin(uDTrec_theta).*cos(bT_phi - uDTrec_phi);
                        P2cosbeta = (3*cosbeta.^2 - 1)/2;

                        Deff = uDTrec_iso.*(1 + 2*bT_Delta.*uDTrec_Delta.*P2cosbeta);
                        Kernel = exp(-bT_trace.*Deff).*exp(-bT_TR1w.*uDTrec_R1)...
                            .*(1 - exp(-bT_TR.*uDTrec_R1)).*exp(-bT_TE.*uDTrec_R2);

                        S_vector_temp = S_vector(indx_BS,1);
                        Kernel_temp = Kernel(indx_BS,:);

    %                     [S_vector_temp,indx] = sort(S_vector_temp,1,'descend');
    %                     Kernel_temp = Kernel_temp(indx,:);

                        uDTrec_w_vtemp = S0*lsqnonneg(Kernel_temp,S_vector_temp/S0);
                        Scalc = Kernel_temp*uDTrec_w_vtemp;
                        chisq = sum((S_vector_temp-Scalc).^2,1)/Npoints;

    %                     figure(2), clf
    %                     plot(1:NpointsBS,Scalc,'k-',1:NpointsBS,S_vector_temp,'b-',[0; NpointsBS],sqrt(chisq)*[1; 1],'k--')
    %                     title([num2str(niter) ' ' num2str(nfail) ' ' num2str(nnonzero) ' ' num2str(chisq,3) ' ' num2str(chisq0,3)])
    %                     pause(.1)

                        niter = niter+1;

                        if chisq > chisq0*(1 - 1e-3)
                            nfail = nfail+1;
                        end

                        if chisq < chisq0*(1 - 0e-2)
                            uDTrec_iso_v = uDTrec_iso_vtemp;
                            uDTrec_Delta_v = uDTrec_Delta_vtemp;
                            uDTrec_theta_v = uDTrec_theta_vtemp;
                            uDTrec_phi_v = uDTrec_phi_vtemp;
                            uDTrec_R1_v = uDTrec_R1_vtemp;
                            uDTrec_R2_v = uDTrec_R2_vtemp;
                            uDTrec_w_v = uDTrec_w_vtemp;
                            chisq0 = chisq;
                        end

                    end


                [uDTrec_w_v,indx] = sort(uDTrec_w_v,'descend');
                uDTrec_iso_v = uDTrec_iso_v(indx);
                uDTrec_Delta_v = uDTrec_Delta_v(indx);
                uDTrec_theta_v = uDTrec_theta_v(indx);
                uDTrec_phi_v = uDTrec_phi_v(indx); 
                uDTrec_R1_v = uDTrec_R1_v(indx); 
                uDTrec_R2_v = uDTrec_R2_v(indx); 
                nnonzero = sum(uDTrec_w_v>0);

                uDTrec_iso_a(:,nBS) = uDTrec_iso_v;
                uDTrec_Delta_a(:,nBS) = uDTrec_Delta_v;
                uDTrec_theta_a(:,nBS) = uDTrec_theta_v;
                uDTrec_phi_a(:,nBS) = uDTrec_phi_v;
                uDTrec_R1_a(:,nBS) = uDTrec_R1_v;
                uDTrec_R2_a(:,nBS) = uDTrec_R2_v;
                uDTrec_w_a(:,nBS) = uDTrec_w_v;

                chisq_v(nBS,1) = chisq;
                niter_v(nBS,1) = niter;
                nfail_v(nBS,1) = nfail;
                nnonzero_v(nBS,1) = nnonzero;
            %     figure(2), clf
            %     plot(1:Npoints,S_vector_temp,'o',1:Npoints,Scalc,'-')
            %     pause(.1)

                p.progress; %Counter for progress report

            end
            p.stop;  
            %toc
            %%

            figure(1), clf
            subplot(2,2,1)
            plot(1:NBS,sqrt(chisq_v),'o')
            ylabel('rms chi'), xlabel('nBS')
            subplot(2,2,2)
            hist(sqrt(chisq_v))
            subplot(2,2,3)
            hist(niter_v)
            subplot(2,2,4)
            hist(nnonzero_v)
            xlabel('nnonzero')

            uDTrec_a.w = uDTrec_w_a;
            uDTrec_a.iso = uDTrec_iso_a;
            uDTrec_a.Delta = uDTrec_Delta_a;
            uDTrec_a.theta = uDTrec_theta_a;
            uDTrec_a.phi = uDTrec_phi_a;
            uDTrec_a.R1 = uDTrec_R1_a;
            uDTrec_a.R2 = uDTrec_R2_a;

            ReconDat = struct('NBS',NBS,'Nnodes',Nnodes,'Dmin',Dmin,'Dmax',Dmax,...
                'R1min',R1min,'R1max',R1max,'R2min',R2min,'R2max',R2max,...
                'nnonzero',nnonzero_v,'chisq',chisq_v);

            save([num2str(expno(nexp)) '/DTR1R2spec'],'uDTrec_a','ReconDat','bT','Signal')
        end
    end
    
end
