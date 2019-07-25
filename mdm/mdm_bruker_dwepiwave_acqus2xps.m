function xps = mdm_bruker_dwepiwave_acqus2xps(data_path, xps_fn)
% function xps = mdm_bruker_dwepiwave_acqus2xps(data_path, xps_fn)
%
% Calculation of b-tensors for gradient waveforms 
% from Topgaard, Phys. Chem. Chem. Phys., 2016.
% http://dx.doi.org/10.1039/c5cp07251d
%
% data_path: directory for gradient text files
% xps_fn (optional) filename of xps

if (nargin < 2), xps_fn = []; end

% Read PV parameters
data_path_pv = [data_path '/'];

if ~strcmp(ReadPVParam(data_path_pv, 'PULPROG'),lower('<rFOV_DWEpiWavev1_04.ppg>')), return, end

DwGradDur = ReadPVParam(data_path_pv, 'DwGradDur')*1e-3 ; 
DwGradTsep = ReadPVParam(data_path_pv, 'DwGradTsep')*1e-3 ;     
DwDelay1 = ReadPVParam(data_path_pv, 'DwDelay1')*1e-3 ;   
DwDelay2 = ReadPVParam(data_path_pv, 'DwDelay2')*1e-3 ;   
DwNB0 = ReadPVParam(data_path_pv, 'DwNB0') ;
DwNDirs = ReadPVParam(data_path_pv, 'DwNDirs') ;
DwNAmplitudes = ReadPVParam(data_path_pv, 'DwNAmplitudes') ;
DwNShapes = ReadPVParam(data_path_pv, 'DwNShapes') ;
DwDir = ReadPVParam(data_path_pv, 'DwDir') ;
DwGradAmpG = ReadPVParam(data_path_pv, 'DwGradAmpG') ;
DwLoopOrder = ReadPVParam(data_path_pv, 'DwLoopOrder') ;
DwGradShapeArray = ReadPVParam(data_path_pv, 'DwGradShapeArray') ;
DwGradShapeStrArr = ReadPVParam(data_path_pv, 'DwGradShapeStrArr') ;
PVM_GradCalConst = ReadPVParam(data_path_pv, 'PVM_GradCalConst') ;
EchoTime=ReadPVParam(data_path_pv, 'EchoTime')*1e-3 ;
RefSliceGrad = ReadPVParam(data_path_pv, 'RefSliceGrad') ;
%RefPul = ReadPVParam(data_path_pv, 'RefPul') ;
SliceSpoilerDur = ReadPVParam(data_path_pv, 'SliceSpoilerDur')*1e-3 ;
SliceSpoilerAmp = ReadPVParam(data_path_pv, 'SliceSpoilerAmp') ;
SliceSpoilerDir = ReadPVParam(data_path_pv, 'SliceSpoilerDir') ; %Seems not to be used. Scope shows DW dirs rotate with slice. 
ACQ_grad_matrix = ReadPVParam(data_path_pv, 'ACQ_grad_matrix') ;    

for nShape = 1:DwNShapes
    parnamR = ['DWGS' num2str(nShape-1) 'R'];
    parvalR = ReadPVParam(data_path_pv, parnamR)' ;
    eval([parnamR ' = parvalR;'])
    parnamP = ['DWGS' num2str(nShape-1) 'P'];
    parvalP = ReadPVParam(data_path_pv, parnamP)' ;
    eval([parnamP ' = parvalP;'])
    parnamS = ['DWGS' num2str(nShape-1) 'S'];
    parvalS = ReadPVParam(data_path_pv, parnamS)' ;
    eval([parnamS ' = parvalS;'])
    %figure(1), clf, plot(1:numel(parvalR),parvalR,'r-', 1:numel(parvalP),parvalP,'g-', 1:numel(parvalS),parvalS,'b-'), pause(1)
end

% Read gradient ramps
% ramp.xa etc maps gradient modulations (a,b,c) to channels (x,y,z)
% Normalized from -1 to +1
ramp.xa = ReadPVParam(data_path_pv, 'DwGAmpRot00')' ;
ramp.xb = ReadPVParam(data_path_pv, 'DwGAmpRot01')' ;
ramp.xc = ReadPVParam(data_path_pv, 'DwGAmpRot02')' ;
ramp.ya = ReadPVParam(data_path_pv, 'DwGAmpRot10')' ;
ramp.yb = ReadPVParam(data_path_pv, 'DwGAmpRot11')' ;
ramp.yc = ReadPVParam(data_path_pv, 'DwGAmpRot12')' ;
ramp.za = ReadPVParam(data_path_pv, 'DwGAmpRot20')' ;
ramp.zb = ReadPVParam(data_path_pv, 'DwGAmpRot21')' ;
ramp.zc = ReadPVParam(data_path_pv, 'DwGAmpRot22')' ;

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
xps = mdm_xps_calc_btpars(xps);

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
    save(fullfile(fileparts(xps_fn),'shapes'),'shapes');    
end
