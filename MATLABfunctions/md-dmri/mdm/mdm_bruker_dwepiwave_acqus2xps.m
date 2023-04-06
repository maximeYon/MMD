function res = mdm_bruker_dwepiwave_acqus2xps(data_path, xps_fn)
% function res = mdm_bruker_dwepiwave_acqus2xps(data_path, xps_fn)
%
% Calculation of b-tensors for gradient waveforms 
% from Matt Budde's sequences rFOV_DWEpiWavev1_04.ppg and mcw_DWEpiWavev7.ppg
%
% data_path: directory for acquisition parameter files
% xps_fn (optional) filename of xps
%

if (nargin < 2), xps_fn = []; end

% Read PV parameters

if ~any(strcmp(mdm_bruker_readpvparam(data_path, 'PULPROG'),{lower('<rFOV_DWEpiWavev1_04.ppg>'); lower('<mcw_DWEpiWavev7.ppg>')})), res = 0; return, end
if ~strcmp(mdm_bruker_readpvparam(data_path, 'ACQ_pipe_status'),lower('Wrapup')), res = 0; return, end


if any(strcmp(mdm_bruker_readpvparam(data_path, 'PULPROG'),{lower('<rFOV_DWEpiWavev1_04.ppg>')}))

    DwGradDur = mdm_bruker_readpvparam(data_path, 'DwGradDur')*1e-3 ; 
    DwGradTsep = mdm_bruker_readpvparam(data_path, 'DwGradTsep')*1e-3 ;     
    DwDelay1 = mdm_bruker_readpvparam(data_path, 'DwDelay1')*1e-3 ;   
    DwDelay2 = mdm_bruker_readpvparam(data_path, 'DwDelay2')*1e-3 ;   
    DwNB0 = mdm_bruker_readpvparam(data_path, 'DwNB0') ;
    DwNDirs = mdm_bruker_readpvparam(data_path, 'DwNDirs') ;
    DwNAmplitudes = mdm_bruker_readpvparam(data_path, 'DwNAmplitudes') ;
    DwNShapes = mdm_bruker_readpvparam(data_path, 'DwNShapes') ;
    DwDir = mdm_bruker_readpvparam(data_path, 'DwDir') ;
    DwGradAmpG = mdm_bruker_readpvparam(data_path, 'DwGradAmpG') ;
    DwLoopOrder = mdm_bruker_readpvparam(data_path, 'DwLoopOrder') ;
    DwGradShapeArray = mdm_bruker_readpvparam(data_path, 'DwGradShapeArray') ;
    DwGradShapeStrArr = mdm_bruker_readpvparam(data_path, 'DwGradShapeStrArr') ;
    PVM_GradCalConst = mdm_bruker_readpvparam(data_path, 'PVM_GradCalConst') ;
    EchoTime=mdm_bruker_readpvparam(data_path, 'EchoTime')*1e-3 ;
    PVM_RepetitionTime=mdm_bruker_readpvparam(data_path, 'PVM_RepetitionTime')*1e-3 ;
    RefSliceGrad = mdm_bruker_readpvparam(data_path, 'RefSliceGrad') ;
    %RefPul = mdm_bruker_readpvparam(data_path, 'RefPul') ;
    SliceSpoilerDur = mdm_bruker_readpvparam(data_path, 'SliceSpoilerDur')*1e-3 ;
    SliceSpoilerAmp = mdm_bruker_readpvparam(data_path, 'SliceSpoilerAmp') ;
    SliceSpoilerDir = mdm_bruker_readpvparam(data_path, 'SliceSpoilerDir') ; %Seems not to be used. Scope shows DW dirs rotate with slice. 
    ACQ_grad_matrix = mdm_bruker_readpvparam(data_path, 'ACQ_grad_matrix') ;    

    for nShape = 1:DwNShapes
        parnamR = ['DWGS' num2str(nShape-1) 'R'];
        parvalR = mdm_bruker_readpvparam(data_path, parnamR)' ;
        eval([parnamR ' = parvalR;'])
        parnamP = ['DWGS' num2str(nShape-1) 'P'];
        parvalP = mdm_bruker_readpvparam(data_path, parnamP)' ;
        eval([parnamP ' = parvalP;'])
        parnamS = ['DWGS' num2str(nShape-1) 'S'];
        parvalS = mdm_bruker_readpvparam(data_path, parnamS)' ;
        eval([parnamS ' = parvalS;'])
        %figure(1), clf, plot(1:numel(parvalR),parvalR,'r-', 1:numel(parvalP),parvalP,'g-', 1:numel(parvalS),parvalS,'b-'), pause(1)
    end

    % Read gradient ramps
    % ramp.xa etc maps gradient modulations (a,b,c) to channels (x,y,z)
    % Normalized from -1 to +1
    ramp.xa = mdm_bruker_readpvparam(data_path, 'DwGAmpRot00')' ;
    ramp.xb = mdm_bruker_readpvparam(data_path, 'DwGAmpRot01')' ;
    ramp.xc = mdm_bruker_readpvparam(data_path, 'DwGAmpRot02')' ;
    ramp.ya = mdm_bruker_readpvparam(data_path, 'DwGAmpRot10')' ;
    ramp.yb = mdm_bruker_readpvparam(data_path, 'DwGAmpRot11')' ;
    ramp.yc = mdm_bruker_readpvparam(data_path, 'DwGAmpRot12')' ;
    ramp.za = mdm_bruker_readpvparam(data_path, 'DwGAmpRot20')' ;
    ramp.zb = mdm_bruker_readpvparam(data_path, 'DwGAmpRot21')' ;
    ramp.zc = mdm_bruker_readpvparam(data_path, 'DwGAmpRot22')' ;

    td1 = numel(ramp.xa);

    % Convert ramps r.ax to uvec axis of q-vector cone
    gpas.avecnorm = sqrt(ramp.xa.^2+ramp.ya.^2+ramp.za.^2);
    gpas.bvecnorm = sqrt(ramp.xb.^2+ramp.yb.^2+ramp.zb.^2);
    gpas.cvecnorm = sqrt(ramp.xc.^2+ramp.yc.^2+ramp.zc.^2);
    gpas.avec = [ramp.xa ramp.ya ramp.za]./repmat(gpas.avecnorm,[1 3]);
    gpas.bvec = [ramp.xb ramp.yb ramp.zb]./repmat(gpas.bvecnorm,[1 3]);
    gpas.cvec = [ramp.xc ramp.yc ramp.zc]./repmat(gpas.cvecnorm,[1 3]);
    gpas.abcrossvec = msf_crossprod_nx3vectors(gpas.avec,gpas.bvec);
    uvec = zeros(td1,3);
    ind = find(gpas.avecnorm>0);
    uvec(ind,:) = gpas.abcrossvec(ind,:);
    ind = find(gpas.cvecnorm>0);
    uvec(ind,:) = gpas.cvec(ind,:);

    gamma = 26.75e7;
    gnams = {'xa','xb','xc','ya','yb','yc','za','zb','zc'};
    for ngnam = 1:numel(gnams)
    gnam = gnams{ngnam};
    eval(['G.' gnam  '= ramp.' gnam '/100*PVM_GradCalConst*2*pi/gamma*1e3;'])    
    end

    Gslice = RefSliceGrad/100*PVM_GradCalConst*2*pi/gamma*1e3;
    Gslicespoil = SliceSpoilerAmp/100*PVM_GradCalConst*2*pi/gamma*1e3;
    slicevector = ACQ_grad_matrix*[0; 0; 1]; % Only an educated guess.

    td1pershape = td1/DwNShapes;
    xps.n = td1;
    xps.b = zeros(td1,1);
    xps.bt = zeros(td1,6);
    xps.u = uvec;
    xps.theta = acos(uvec(:,3));
    xps.phi = atan2(uvec(:,2),uvec(:,1));
    shapes = cell(DwNShapes,1);
    for nShape = 1:DwNShapes
        if strcmp(DwLoopOrder,'dir_amp_shape')
            ind = (1:(td1pershape)) + (nShape-1)*td1pershape;
        else
            warning(['Code currently only works for DwLoopOrder = dir_amp_shape.']);
        end

        eval(['G.a = DWGS' num2str(nShape-1) 'R;'])
        eval(['G.b = DWGS' num2str(nShape-1) 'P;'])
        eval(['G.c = DWGS' num2str(nShape-1) 'S;'])
        shapes{nShape,1} = [G.a G.b G.c];

        Ndt = numel(G.a);
        dt = DwGradDur/Ndt;

        Ndt_GradTsep = round(DwGradTsep/dt); Ndt_GradTsep = 2*ceil(Ndt_GradTsep/2);
        Ndt_SliceSpoiler = round(SliceSpoilerDur/dt);
        Ndt_Slice = round((DwGradTsep-2*SliceSpoilerDur-DwDelay1-DwDelay2)/dt); Ndt_Slice = 2*ceil(Ndt_Slice/2);
        Ndt_G0 = round((Ndt_GradTsep - 2*Ndt_SliceSpoiler - Ndt_Slice)/2);

        ind_slice = logical([zeros(Ndt,1); ones(Ndt_GradTsep,1); zeros(Ndt,1)]);
        Ndt = numel(ind_slice);
        ind_pre180 = logical([ones(Ndt/2,1); zeros(Ndt/2,1)]); 
        ind_pre180 = repmat(ind_pre180,[1, td1pershape]);

        G.slice = [zeros(Ndt_G0,1); Gslicespoil*ones(Ndt_SliceSpoiler,1); Gslice*ones(Ndt_Slice,1); ...
            Gslicespoil*ones(Ndt_SliceSpoiler,1); zeros(Ndt_G0,1)];        

        G.a = [G.a; zeros(Ndt_GradTsep,1); G.a];
        G.b = [G.b; zeros(Ndt_GradTsep,1); G.b];
        G.c = [G.c; zeros(Ndt_GradTsep,1); G.c];

        %figure(1), clf, plot(dt*(1:Ndt),G.a,'r-',dt*(1:Ndt),G.b,'g-',dt*(1:Ndt),G.c,'b-'), pause(1), return       
        if Ndt > 1
            G.x = repmat(G.a,[1, td1pershape]).*repmat(G.xa(ind)',[Ndt, 1]) + ...
                repmat(G.b,[1, td1pershape]).*repmat(G.xb(ind)',[Ndt, 1]) + ...
                repmat(G.c,[1, td1pershape]).*repmat(G.xc(ind)',[Ndt, 1]);
            G.y = repmat(G.a,[1, td1pershape]).*repmat(G.ya(ind)',[Ndt, 1]) + ...
                repmat(G.b,[1, td1pershape]).*repmat(G.yb(ind)',[Ndt, 1]) + ...
                repmat(G.c,[1, td1pershape]).*repmat(G.yc(ind)',[Ndt, 1]);
            G.z = repmat(G.a,[1, td1pershape]).*repmat(G.za(ind)',[Ndt, 1]) + ...
                repmat(G.b,[1, td1pershape]).*repmat(G.zb(ind)',[Ndt, 1]) + ...
                repmat(G.c,[1, td1pershape]).*repmat(G.zc(ind)',[Ndt, 1]);

            G.x(ind_slice,:) = slicevector(1)*repmat(G.slice,[1, td1pershape]);
            G.y(ind_slice,:) = slicevector(2)*repmat(G.slice,[1, td1pershape]);
            G.z(ind_slice,:) = slicevector(3)*repmat(G.slice,[1, td1pershape]);

            G.x(ind_pre180) = -G.x(ind_pre180);
            G.y(ind_pre180) = -G.y(ind_pre180);
            G.z(ind_pre180) = -G.z(ind_pre180);

            %figure(1), clf, plot(dt*(1:Ndt)',G.x,'r-',dt*(1:Ndt)',G.y,'g-',dt*(1:Ndt)',G.z,'b-'), pause(.1)        

            % Dephasing vector F.x in SI
            F.x = cumsum(G.x*dt,1);
            F.y = cumsum(G.y*dt,1);
            F.z = cumsum(G.z*dt,1);
            F.r = sqrt(F.x.^2 + F.y.^2 + F.z.^2); % Magnitude

            % Diffusion weighting matrix b                
            bt.xx = gamma^2*sum(F.x.*F.x*dt,1)';
            bt.xy = gamma^2*sum(F.x.*F.y*dt,1)';
            bt.xz = gamma^2*sum(F.x.*F.z*dt,1)';
            bt.yy = gamma^2*sum(F.y.*F.y*dt,1)';
            bt.yz = gamma^2*sum(F.y.*F.z*dt,1)';
            bt.zz = gamma^2*sum(F.z.*F.z*dt,1)';

            b = (bt.xx + bt.yy + bt.zz); % trace

            % Save as xps
            xps.b(ind,1) = b;
            xps.bt(ind,:) = [bt.xx bt.yy bt.zz sqrt(2)*[bt.xy bt.xz bt.yz]];
            %xps.u(ind,:) = [zeros(DwNB0,3); repmat(DwDir,[DwNAmplitudes 1])];
        end
    end
    
    xps.te = EchoTime*ones(td1,1);
    xps.tr = PVM_RepetitionTime*ones(td1,1);
    xps_temp = mdm_xps_calc_btpars(xps);
    xps.b_delta = xps_temp.b_delta;
    xps.b_eta = xps_temp.b_eta;

    % figure(2), clf
    % subplot(2,1,1)
    % hx = plot(1:xps.n,G.xa,'ro',1:xps.n,G.xb,'rs',1:xps.n,G.xc,'rx');
    % hold on
    % hy = plot(1:xps.n,G.ya,'go',1:xps.n,G.yb,'gs',1:xps.n,G.yc,'gx');
    % hz = plot(1:xps.n,G.za,'bo',1:xps.n,G.zb,'bs',1:xps.n,G.zc,'bx');
    % set(hx,'MarkerSize',10), set(hy,'MarkerSize',8), set(hz,'MarkerSize',6)
    % title(['td1 = ' num2str(xps.n)])
    % subplot(2,1,2)
    % ha = plot(1:xps.n,sqrt(G.xa.^2+G.ya.^2+G.za.^2),'ro');
    % hold on
    % hb = plot(1:xps.n,sqrt(G.xb.^2+G.yb.^2+G.zb.^2),'gs');
    % hc = plot(1:xps.n,sqrt(G.xc.^2+G.yc.^2+G.zc.^2),'bx');

    if (~isempty(xps_fn))
        msf_mkdir(fileparts(xps_fn));
        save(xps_fn,'xps');
        %save(fullfile(fileparts(xps_fn),'shapes'),'shapes');    
    end

    res = 1;
    
elseif any(strcmp(mdm_bruker_readpvparam(data_path, 'PULPROG'),{lower('<mcw_DWEpiWavev7.ppg>')}))
    DwGradDur1 = mdm_bruker_readpvparam(data_path, 'DwGradDur1')*1e-3 ; 
    DwGradDur2 = mdm_bruker_readpvparam(data_path, 'DwGradDur2')*1e-3 ; 
    DwGradTsep = mdm_bruker_readpvparam(data_path, 'DwGradTsep')*1e-3 ;     
    DwDelay1 = mdm_bruker_readpvparam(data_path, 'DwDelay1')*1e-3 ;   
    DwDelay2 = mdm_bruker_readpvparam(data_path, 'DwDelay2')*1e-3 ;   
    DwNB0 = mdm_bruker_readpvparam(data_path, 'DwNB0') ;
    DwNDirs = mdm_bruker_readpvparam(data_path, 'DwNDirs') ;
    DwNAmplitudes = mdm_bruker_readpvparam(data_path, 'DwNAmplitudes') ;
    DwNShapes = mdm_bruker_readpvparam(data_path, 'DwNShapes') ;
    DwDir = mdm_bruker_readpvparam(data_path, 'DwDir') ;
    DwGradAmpG = mdm_bruker_readpvparam(data_path, 'DwGradAmpG') ;
    DwLoopOrder = mdm_bruker_readpvparam(data_path, 'DwLoopOrder') ;
    DwGradShapeArray1 = mdm_bruker_readpvparam(data_path, 'DwGradShapeArray1') ;
    DwGradShapeArray2 = mdm_bruker_readpvparam(data_path, 'DwGradShapeArray2') ;
    DwGradShapeStrArr1 = mdm_bruker_readpvparam(data_path, 'DwGradShapeStrArr1') ;
    DwGradShapeStrArr1 = mdm_bruker_readpvparam(data_path, 'DwGradShapeStrArr1') ;
    PVM_GradCalConst = mdm_bruker_readpvparam(data_path, 'PVM_GradCalConst') ;
    EchoTime=mdm_bruker_readpvparam(data_path, 'EchoTime')*1e-3 ;
    PVM_RepetitionTime=mdm_bruker_readpvparam(data_path, 'PVM_RepetitionTime')*1e-3 ;
    RefSliceGrad = mdm_bruker_readpvparam(data_path, 'RefSliceGrad') ;
    %RefPul = mdm_bruker_readpvparam(data_path, 'RefPul') ;
    SliceSpoilerDur = mdm_bruker_readpvparam(data_path, 'SliceSpoilerDur')*1e-3 ;
    SliceSpoilerAmp = mdm_bruker_readpvparam(data_path, 'SliceSpoilerAmp') ;
    SliceSpoilerDir = mdm_bruker_readpvparam(data_path, 'SliceSpoilerDir') ; %Seems not to be used. Scope shows DW dirs rotate with slice. 
    ACQ_grad_matrix = mdm_bruker_readpvparam(data_path, 'ACQ_grad_matrix') ;    

    for nShape = 1:DwNShapes
        parnamR = ['DW1GS' num2str(nShape-1) 'R'];
        parvalR1 = mdm_bruker_readpvparam(data_path, parnamR)' ;
        eval([parnamR ' = parvalR1;'])
        parnamP = ['DW1GS' num2str(nShape-1) 'P'];
        parvalP1 = mdm_bruker_readpvparam(data_path, parnamP)' ;
        eval([parnamP ' = parvalP1;'])
        parnamS = ['DW1GS' num2str(nShape-1) 'S'];
        parvalS1 = mdm_bruker_readpvparam(data_path, parnamS)' ;
        eval([parnamS ' = parvalS1;'])
        parnamR = ['DW2GS' num2str(nShape-1) 'R'];
        parvalR2 = mdm_bruker_readpvparam(data_path, parnamR)' ;
        eval([parnamR ' = parvalR2;'])
        parnamP = ['DW2GS' num2str(nShape-1) 'P'];
        parvalP2 = mdm_bruker_readpvparam(data_path, parnamP)' ;
        eval([parnamP ' = parvalP2;'])
        parnamS = ['DW2GS' num2str(nShape-1) 'S'];
        parvalS2 = mdm_bruker_readpvparam(data_path, parnamS)' ;
        eval([parnamS ' = parvalS2;'])
        %figure(1), clf, subplot(1,2,1), plot(1:numel(parvalR1),parvalR1,'r-', 1:numel(parvalP1),parvalP1,'g-', 1:numel(parvalS1),parvalS1,'b-'), subplot(1,2,2), plot(1:numel(parvalR2),parvalR2,'r-', 1:numel(parvalP2),parvalP2,'g-', 1:numel(parvalS2),parvalS2,'b-'), pause(1)
    end

    % Read gradient ramps
    % ramp.xa etc maps gradient modulations (a,b,c) to channels (x,y,z)
    % Normalized from -1 to +1
    ramp.xa = mdm_bruker_readpvparam(data_path, 'DwGAmpRot00')' ;
    ramp.xb = mdm_bruker_readpvparam(data_path, 'DwGAmpRot01')' ;
    ramp.xc = mdm_bruker_readpvparam(data_path, 'DwGAmpRot02')' ;
    ramp.ya = mdm_bruker_readpvparam(data_path, 'DwGAmpRot10')' ;
    ramp.yb = mdm_bruker_readpvparam(data_path, 'DwGAmpRot11')' ;
    ramp.yc = mdm_bruker_readpvparam(data_path, 'DwGAmpRot12')' ;
    ramp.za = mdm_bruker_readpvparam(data_path, 'DwGAmpRot20')' ;
    ramp.zb = mdm_bruker_readpvparam(data_path, 'DwGAmpRot21')' ;
    ramp.zc = mdm_bruker_readpvparam(data_path, 'DwGAmpRot22')' ;

    td1 = numel(ramp.xa);

    % Convert ramps r.ax to uvec axis of q-vector cone
    gpas.avecnorm = sqrt(ramp.xa.^2+ramp.ya.^2+ramp.za.^2);
    gpas.bvecnorm = sqrt(ramp.xb.^2+ramp.yb.^2+ramp.zb.^2);
    gpas.cvecnorm = sqrt(ramp.xc.^2+ramp.yc.^2+ramp.zc.^2);
    gpas.avec = [ramp.xa ramp.ya ramp.za]./repmat(gpas.avecnorm,[1 3]);
    gpas.bvec = [ramp.xb ramp.yb ramp.zb]./repmat(gpas.bvecnorm,[1 3]);
    gpas.cvec = [ramp.xc ramp.yc ramp.zc]./repmat(gpas.cvecnorm,[1 3]);
    gpas.abcrossvec = msf_crossprod_nx3vectors(gpas.avec,gpas.bvec);
    uvec = zeros(td1,3);
    ind = find(gpas.avecnorm>0);
    uvec(ind,:) = gpas.abcrossvec(ind,:);
    ind = find(gpas.cvecnorm>0);
    uvec(ind,:) = gpas.cvec(ind,:);

    gamma = 26.75e7;
    gnams = {'xa','xb','xc','ya','yb','yc','za','zb','zc'};
    for ngnam = 1:numel(gnams)
    gnam = gnams{ngnam};
    eval(['G.' gnam  '= ramp.' gnam '/100*PVM_GradCalConst*2*pi/gamma*1e3;'])    
    end

    Gslice = RefSliceGrad/100*PVM_GradCalConst*2*pi/gamma*1e3;
    Gslicespoil = SliceSpoilerAmp/100*PVM_GradCalConst*2*pi/gamma*1e3;
    slicevector = ACQ_grad_matrix*[0; 0; 1]; % Only an educated guess.

    td1pershape = td1/DwNShapes;
    xps.n = td1;
    xps.b = zeros(td1,1);
    xps.bt = zeros(td1,6);
    xps.u = uvec;
    xps.theta = acos(uvec(:,3));
    xps.phi = atan2(uvec(:,2),uvec(:,1));
    shapes1 = cell(DwNShapes,1);
    shapes2 = cell(DwNShapes,1);
    for nShape = 1:DwNShapes
        if strcmp(DwLoopOrder,'dir_amp_shape')
            ind = (1:(td1pershape)) + (nShape-1)*td1pershape;
        else
            warning(['Code currently only works for DwLoopOrder = dir_amp_shape.']);
        end

        eval(['G1.a = DW1GS' num2str(nShape-1) 'R;'])
        eval(['G1.b = DW1GS' num2str(nShape-1) 'P;'])
        eval(['G1.c = DW1GS' num2str(nShape-1) 'S;'])
        shapes1{nShape,1} = [G1.a G1.b G1.c];

        Ndt1 = numel(G1.a);
        dt1 = DwGradDur1/Ndt1;

        eval(['G2.a = DW2GS' num2str(nShape-1) 'R;'])
        eval(['G2.b = DW2GS' num2str(nShape-1) 'P;'])
        eval(['G2.c = DW2GS' num2str(nShape-1) 'S;'])
        shapes2{nShape,1} = [G2.a G2.b G2.c];

        Ndt2 = numel(G2.a);
        dt2 = DwGradDur2/Ndt2;

        Ndt_GradTsep = round(DwGradTsep/dt1); Ndt_GradTsep = 2*ceil(Ndt_GradTsep/2);
        Ndt_SliceSpoiler = round(SliceSpoilerDur/dt1);
        Ndt_Slice = round((DwGradTsep-2*SliceSpoilerDur-DwDelay1-DwDelay2)/dt1); Ndt_Slice = 2*ceil(Ndt_Slice/2);
        Ndt_G0 = round((Ndt_GradTsep - 2*Ndt_SliceSpoiler - Ndt_Slice)/2);

        ind_slice = logical([zeros(Ndt1,1); ones(Ndt_GradTsep,1); zeros(Ndt2,1)]);
        dt_vector = [dt1*ones(Ndt1,1); dt1*ones(Ndt_GradTsep,1); dt2*ones(Ndt2,1)];
        
        Ndt = numel(ind_slice);
        ind_pre180 = logical([ones(Ndt1,1); ones(Ndt_GradTsep/2,1); zeros(Ndt_GradTsep/2,1); zeros(Ndt2,1)]);
        ind_pre180 = repmat(ind_pre180,[1, td1pershape]);

        G.slice = [zeros(Ndt_G0,1); Gslicespoil*ones(Ndt_SliceSpoiler,1); Gslice*ones(Ndt_Slice,1); ...
            Gslicespoil*ones(Ndt_SliceSpoiler,1); zeros(Ndt_G0,1)];        

        G.a = [G1.a; zeros(Ndt_GradTsep,1); G2.a];
        G.b = [G1.b; zeros(Ndt_GradTsep,1); G2.b];
        G.c = [G1.c; zeros(Ndt_GradTsep,1); G2.c];

        %figure(1), clf, plot(dt*(1:Ndt),G.a,'r-',dt*(1:Ndt),G.b,'g-',dt*(1:Ndt),G.c,'b-'), pause(1), return       
        if Ndt > 1
            G.x = repmat(G.a,[1, td1pershape]).*repmat(G.xa(ind)',[Ndt, 1]) + ...
                repmat(G.b,[1, td1pershape]).*repmat(G.xb(ind)',[Ndt, 1]) + ...
                repmat(G.c,[1, td1pershape]).*repmat(G.xc(ind)',[Ndt, 1]);
            G.y = repmat(G.a,[1, td1pershape]).*repmat(G.ya(ind)',[Ndt, 1]) + ...
                repmat(G.b,[1, td1pershape]).*repmat(G.yb(ind)',[Ndt, 1]) + ...
                repmat(G.c,[1, td1pershape]).*repmat(G.yc(ind)',[Ndt, 1]);
            G.z = repmat(G.a,[1, td1pershape]).*repmat(G.za(ind)',[Ndt, 1]) + ...
                repmat(G.b,[1, td1pershape]).*repmat(G.zb(ind)',[Ndt, 1]) + ...
                repmat(G.c,[1, td1pershape]).*repmat(G.zc(ind)',[Ndt, 1]);

            G.x(ind_slice,:) = slicevector(1)*repmat(G.slice,[1, td1pershape]);
            G.y(ind_slice,:) = slicevector(2)*repmat(G.slice,[1, td1pershape]);
            G.z(ind_slice,:) = slicevector(3)*repmat(G.slice,[1, td1pershape]);

            G.x(ind_pre180) = -G.x(ind_pre180);
            G.y(ind_pre180) = -G.y(ind_pre180);
            G.z(ind_pre180) = -G.z(ind_pre180);

            %figure(1), clf, plot(cumsum(dt_vector,1),G.x,'r-',cumsum(dt_vector,1),G.y,'g-',cumsum(dt_vector,1),G.z,'b-'), pause(.1)        

            % Dephasing vector F.x in SI
            F.x = cumsum(G.x.*dt_vector,1);
            F.y = cumsum(G.y.*dt_vector,1);
            F.z = cumsum(G.z.*dt_vector,1);
            F.r = sqrt(F.x.^2 + F.y.^2 + F.z.^2); % Magnitude

            % Diffusion weighting matrix b                
            bt.xx = gamma^2*sum(F.x.*F.x.*dt_vector,1)';
            bt.xy = gamma^2*sum(F.x.*F.y.*dt_vector,1)';
            bt.xz = gamma^2*sum(F.x.*F.z.*dt_vector,1)';
            bt.yy = gamma^2*sum(F.y.*F.y.*dt_vector,1)';
            bt.yz = gamma^2*sum(F.y.*F.z.*dt_vector,1)';
            bt.zz = gamma^2*sum(F.z.*F.z.*dt_vector,1)';

            b = (bt.xx + bt.yy + bt.zz); % trace

            % Save as xps
            xps.b(ind,1) = b;
            xps.bt(ind,:) = [bt.xx bt.yy bt.zz sqrt(2)*[bt.xy bt.xz bt.yz]];
            %xps.u(ind,:) = [zeros(DwNB0,3); repmat(DwDir,[DwNAmplitudes 1])];
        end
    end

    xps.te = EchoTime*ones(td1,1);
    xps.tr = PVM_RepetitionTime*ones(td1,1);
    xps_temp = mdm_xps_calc_btpars(xps);
    xps.b_delta = xps_temp.b_delta;
    xps.b_eta = xps_temp.b_eta;

%     figure(2), clf
%     subplot(2,1,1)
%     hx = plot(1:xps.n,G.xa,'ro',1:xps.n,G.xb,'rs',1:xps.n,G.xc,'rx');
%     hold on
%     hy = plot(1:xps.n,G.ya,'go',1:xps.n,G.yb,'gs',1:xps.n,G.yc,'gx');
%     hz = plot(1:xps.n,G.za,'bo',1:xps.n,G.zb,'bs',1:xps.n,G.zc,'bx');
%     set(hx,'MarkerSize',10), set(hy,'MarkerSize',8), set(hz,'MarkerSize',6)
%     title(['td1 = ' num2str(xps.n)])
%     subplot(2,1,2)
%     ha = plot(1:xps.n,sqrt(G.xa.^2+G.ya.^2+G.za.^2),'ro');
%     hold on
%     hb = plot(1:xps.n,sqrt(G.xb.^2+G.yb.^2+G.zb.^2),'gs');
%     hc = plot(1:xps.n,sqrt(G.xc.^2+G.yc.^2+G.zc.^2),'bx');

    if (~isempty(xps_fn))
        msf_mkdir(fileparts(xps_fn));
        save(xps_fn,'xps');
        %save(fullfile(fileparts(xps_fn),'shapes'),'shapes');    
    end

    res = 1;
end
