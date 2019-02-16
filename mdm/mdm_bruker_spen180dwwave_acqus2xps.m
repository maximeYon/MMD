function xps = mdm_bruker_spen180dwwave_acqus2xps(data_path, xps_fn)
% function xps = mdm_bruker_dt_axderare2d_acqus2xps(data_path, xps_fn)
%
% Calculation of b-tensors for the Bruker pulse program DT_axderare2d
% as used in Topgaard, Phys. Chem. Chem. Phys., 2016.
% http://dx.doi.org/10.1039/c5cp07251d
%
% data_path: directory for gradient text files
% xps_fn (optional) filename of xps

if (nargin < 2), xps_fn = []; end

data_path_pv = [data_path '/'];

DwGradDur = ReadPVParam(data_path_pv, 'DwGradDur')*1e-3 ;
DwNShapes = ReadPVParam(data_path_pv, 'DwNShapes') ;
DwNDirs = ReadPVParam(data_path_pv, 'DwNDirs') ;
DwNB0 = ReadPVParam(data_path_pv, 'DwNB0') ;
DwDir = ReadPVParam(data_path_pv, 'DwDir') ;
DwGradAmpG = ReadPVParam(data_path_pv, 'DwGradAmpG') ;
DwLoopOrder = ReadPVParam(data_path_pv, 'DwLoopOrder') ;
DwGradShapeArray = ReadPVParam(data_path_pv, 'DwGradShapeArray') ;
PVM_GradCalConst = ReadPVParam(data_path_pv, 'PVM_GradCalConst') ;
DwNDiffExp=ReadPVParam(data_path_pv, 'DwNDiffExp') ;
EchoTime=ReadPVParam(data_path_pv, 'EchoTime') ;
EffectiveEchoTime=ReadPVParam(data_path_pv, 'EffectiveEchoTime')*1e-3 ;
MultipleEchos=ReadPVParam(data_path_pv, 'MultipleEchos') ;

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

% G.xa = ramp.xa/100*PVM_GradCalConst*2*pi/gamma*1e3;
% G.xb = ramp.xb/100*PVM_GradCalConst*2*pi/gamma*1e3;
% G.xc = ramp.xc/100*PVM_GradCalConst*2*pi/gamma*1e3;
% G.ya = ramp.ya/100*PVM_GradCalConst*2*pi/gamma*1e3;
% G.yb = ramp.yb/100*PVM_GradCalConst*2*pi/gamma*1e3;
% G.yc = ramp.yc/100*PVM_GradCalConst*2*pi/gamma*1e3;
% G.za = ramp.za/100*PVM_GradCalConst*2*pi/gamma*1e3;
% G.zb = ramp.zb/100*PVM_GradCalConst*2*pi/gamma*1e3;
% G.zc = ramp.zc/100*PVM_GradCalConst*2*pi/gamma*1e3;

td1pershape = td1/DwNShapes;
xps.n = td1;
xps.b = zeros(td1,1);
xps.bt = zeros(td1,6);
xps.u = uvec;
xps.theta = acos(uvec(:,3));
xps.phi = atan2(uvec(:,2),uvec(:,1));

for nShape = 1:DwNShapes
    ind = (1:(td1pershape)) + (nShape-1)*td1pershape;
    eval(['G.a = DWGS' num2str(nShape-1) 'R;'])
    eval(['G.b = DWGS' num2str(nShape-1) 'P;'])
    eval(['G.c = DWGS' num2str(nShape-1) 'S;'])
    Ndt = numel(G.a);
    dt = DwGradDur/Ndt;
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
    end
end

xps_old = xps;
xps = mdm_xps_calc_btpars(xps_old);


xps.te = EchoTime*ones(xps.n,1);

if (~isempty(DwNDiffExp) && strcmp(MultipleEchos,'yes'))

    Nte = numel(EffectiveEchoTime);
    xps_cell = cell(Nte,1);
    for nte = 1:Nte
        xps_cell{nte} = xps;
        xps_cell{nte}.te =EffectiveEchoTime(nte)*ones(xps.n,1);
    end
    xps = mdm_xps_merge(xps_cell);
    

    EpiAcquistionTime = ReadPVParam(data_path_pv, 'EpiAcquistionTime')*1e-3 ;
    PVM_Matrix = ReadPVParam(data_path_pv, 'PVM_Matrix');

    %xps.dte = repmat(EpiAcquistionTime*fliplr(-1/2:1/(PVM_Matrix(1,2)-1):1/2),[xps.n 1]);
    xps.dte = repmat(EpiAcquistionTime*(-1/2:1/(PVM_Matrix(1,2)-1):1/2),[xps.n 1]);

end

% figure(1), clf
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
% 
% figure(2), clf
% subplot(5,1,1)
% hx = plot(1:xps.n,xps.b,'ro');
% title(['td1 = ' num2str(xps.n)])
% ylabel('b')
% subplot(5,1,2)
% ha = plot(1:xps.n,xps.b_delta,'ro');
% ylabel('b_\Delta')
% subplot(5,1,3)
% ha = plot(1:xps.n,xps.theta,'ro');
% ylabel('\theta')
% subplot(5,1,4)
% ha = plot(1:xps.n,xps.phi,'ro');
% ylabel('\phi')
% subplot(5,1,5)
% ha = plot(1:xps.n,xps.te,'ro');
% ylabel('TE')

if (~isempty(xps_fn)), save(xps_fn,'xps'); end

