function xps = mdm_bruker_diff_wave_sym_msme_acqus2xps(data_path, xps_fn, rps)
% function xps = mdm_bruker_diff_wave_sym_msme_acqus2xps(data_path, xps_fn)
%
% Calculation of b-tensors for gradient waveforms 
% from Topgaard, Phys. Chem. Chem. Phys., 2016.
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
PVM_EchoTime=mdm_bruker_readpvparam(data_path, 'PVM_EchoTime')*1e-3 ;
PVM_RepetitionTime=mdm_bruker_readpvparam(data_path, 'PVM_RepetitionTime')*1e-3 ;
EffectiveEchoTime=ReadPVParam(data_path_pv, 'EffectiveEchoTime')*1e-3 ;
MultipleEchos=ReadPVParam(data_path_pv, 'MultipleEchos') ;

DwGradTsep = mdm_bruker_readpvparam(data_path, 'DwGradTsep')*1e-3;  
gradAmpl=ReadPVParam(data_path_pv, 'ACQ_gradient_amplitude') ;
Delays=ReadPVParam(data_path_pv, 'D') ; % in seconde
Pulses=ReadPVParam(data_path_pv, 'P') ; % in useconde

g0 = gradAmpl(0+1);
p0 = Pulses(0+1)*1e-6;
d2 = Delays(2+1);
g2 = gradAmpl(2+1);
g7 = gradAmpl(7+1);
d1 = Delays(1+1);
DwD1 = ReadPVParam(data_path_pv, 'DwD1')*1e-3; % Not defined?
d3 = Delays(3+1);
p1 = Pulses(1+1)*1e-6;
d5 = Delays(5+1);
g8 = gradAmpl(8+1);
g3 = gradAmpl(3+1);
d6 = Delays(6+1);
d3 = Delays(3+1);
de=ReadPVParam(data_path_pv, 'DE')*1e-6;
denab = d3 - de;
g5 = gradAmpl(5+1);
depa=ReadPVParam(data_path_pv, 'DEPA')*1e-6; % Not defined?
rdepa=de; % Should be rdepa=de-depa

sw_h=ReadPVParam(data_path_pv, 'SW_h');
acq_size=ReadPVParam(data_path_pv, 'ACQ_size');
td = acq_size(1);
dw = .5/sw_h;

aqq=dw*td;

% nShape = 1;
% eval(['G.a = DWGS' num2str(nShape-1) 'R;'])
% Ndt = numel(G.a);
% if Ndt == 1
%     Ndt = round(DwGradDur/8e-6);
% end

dt = .5*8e-6;
Ndt_DwGradDur = round(DwGradDur/dt);

for nShape = 1:DwNShapes
    parnamR = ['DWGS' num2str(nShape-1) 'R'];
    parvalR = ReadPVParam(data_path_pv, parnamR)' ;
    parnamP = ['DWGS' num2str(nShape-1) 'P'];
    parvalP = ReadPVParam(data_path_pv, parnamP)' ;
    parnamS = ['DWGS' num2str(nShape-1) 'S'];
    parvalS = ReadPVParam(data_path_pv, parnamS)' ;
%     figure(1), clf, plot(1:numel(parvalR),parvalR,'r-', 1:numel(parvalP),parvalP,'g-', 1:numel(parvalS),parvalS,'b-'), pause%(1)

    if numel(parvalR) == 1
        Ndt = round(DwGradDur/dt);
        parvalR = parvalR*ones(Ndt,1);
        parvalP = parvalP*ones(Ndt,1);
        parvalS = parvalS*ones(Ndt,1);
    end
        
    parvalR = interp1(linspace(0,1,numel(parvalR))',parvalR,linspace(0,1,Ndt_DwGradDur)','makima');
    parvalP = interp1(linspace(0,1,numel(parvalP))',parvalP,linspace(0,1,Ndt_DwGradDur)','makima');
    parvalS = interp1(linspace(0,1,numel(parvalS))',parvalS,linspace(0,1,Ndt_DwGradDur)','makima');

    eval([parnamR ' = parvalR;'])
    eval([parnamP ' = parvalP;'])
    eval([parnamS ' = parvalS;'])
    
%     figure(1), clf, plot(1:numel(parvalR),parvalR,'r-', 1:numel(parvalP),parvalP,'g-', 1:numel(parvalS),parvalS,'b-'), pause%(1)
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
shapes = cell(DwNShapes,1);

Ndt_p0half = round(p0/2/dt);
Ndt_d2 = round(d2/dt);
Ndt_d1 = round(d1/dt);
Ndt_d3 = round(d3/dt);
Ndt_p1half = round(p1/2/dt);
Ndt_d5 = round(d5/dt);
Ndt_d6 = round(d6/dt);
Ndt_denab = round(denab/dt);
Ndt_rdepa = round(rdepa/dt);
Ndt_aqqhalf = round(aqq/2/dt);


Ndt_GradTsep = round(DwGradTsep/dt); Ndt_GradTsep = 2*ceil(Ndt_GradTsep/2);
% Ndt_tot = 2*Ndt+Ndt_GradTsep;
Ndt_preref = Ndt_p0half + Ndt_d2 + Ndt_d1 + ...
    Ndt_DwGradDur + ...
    Ndt_d3 + Ndt_p1half;
Ndt_postref = Ndt_p1half + Ndt_d5 + Ndt_d6 + ...
    Ndt_DwGradDur + ...
    Ndt_denab + Ndt_rdepa + Ndt_aqqhalf;
Ndt_tot = Ndt_preref + Ndt_postref;

xps.gwf_x = zeros(td1,Ndt_tot); %Assumes all shapes have the same number of time steps
xps.gwf_y = zeros(td1,Ndt_tot);
xps.gwf_z = zeros(td1,Ndt_tot);

unit2si = 1/100*PVM_GradCalConst*2*pi/gamma*1e3;

G.exc_x = [0*ones(Ndt_p0half,1); g2*ones(Ndt_d2,1); 0*ones(Ndt_d1,1)]*unit2si;
G.exc_y = [0*ones(Ndt_p0half,1); 0*ones(Ndt_d2,1); 0*ones(Ndt_d1,1)]*unit2si;
G.exc_z = [g0*ones(Ndt_p0half,1); g7*ones(Ndt_d2,1); 0*ones(Ndt_d1,1)]*unit2si;

G.ref_x = [0*ones(Ndt_d3,1); 0*ones(Ndt_p1half,1); 0*ones(Ndt_p1half,1); g8*ones(Ndt_d5,1); 0*ones(Ndt_d6,1)]*unit2si;
G.ref_y = [0*ones(Ndt_d3,1); 0*ones(Ndt_p1half,1); 0*ones(Ndt_p1half,1); 0*ones(Ndt_d5,1); 0*ones(Ndt_d6,1)]*unit2si;
G.ref_z = [g0*ones(Ndt_d3,1); g0*ones(Ndt_p1half,1); g0*ones(Ndt_p1half,1); g3*ones(Ndt_d5,1); 0*ones(Ndt_d6,1)]*unit2si;

G.acq_x = [g5*ones(Ndt_denab,1); g5*ones(Ndt_rdepa,1); g5*ones(Ndt_aqqhalf,1)]*unit2si;
G.acq_y = [0*ones(Ndt_denab,1); 0*ones(Ndt_rdepa,1); 0*ones(Ndt_aqqhalf,1)]*unit2si;
G.acq_z = [0*ones(Ndt_denab,1); 0*ones(Ndt_rdepa,1); 0*ones(Ndt_aqqhalf,1)]*unit2si;

G.im_x = [G.exc_x; zeros(Ndt_DwGradDur,1); G.ref_x; zeros(Ndt_DwGradDur,1); G.acq_x];
G.im_y = [G.exc_y; zeros(Ndt_DwGradDur,1); G.ref_y; zeros(Ndt_DwGradDur,1); G.acq_y];
G.im_z = [G.exc_z; zeros(Ndt_DwGradDur,1); G.ref_z; zeros(Ndt_DwGradDur,1); G.acq_z];

% figure(1), clf, plot(dt*(1:Ndt_tot)',G.im_x,'r-',dt*(1:Ndt_tot)',G.im_y,'g-',dt*(1:Ndt_tot)',G.im_z,'b-'), %return 

ind_im = logical([ones(numel(G.exc_x),1); zeros(Ndt_DwGradDur,1); ones(numel(G.ref_x),1); zeros(Ndt_DwGradDur,1); ones(numel(G.acq_x),1)]);
ind_preref = logical([ones(Ndt_preref,1); zeros(Ndt_postref,1)]);

fh_shapes = figure(1); clf
axh_grad = axes(fh_shapes,'position',[.1 .73 .85 .25]);
axh_deph = axes(fh_shapes,'position',[.1 .45 .85 .25]);
axh_spec = axes(fh_shapes,'position',[.1 .1 .85 .25]);
hold(axh_grad,'on')
hold(axh_deph,'on')
hold(axh_spec,'on')
for nShape = 1:DwNShapes
    ind = (1:(td1pershape)) + (nShape-1)*td1pershape;
    eval(['G.a = DWGS' num2str(nShape-1) 'R;'])
    eval(['G.b = DWGS' num2str(nShape-1) 'P;'])
    eval(['G.c = DWGS' num2str(nShape-1) 'S;'])
    
    Ndt = numel(G.a);
%     if Ndt == 1
%         Ndt = round(DwGradDur/8e-6);
%         G.a = zeros(Ndt,1);
%         G.b = zeros(Ndt,1);
%         G.c = zeros(Ndt,1);
%     end
%     dt = DwGradDur/Ndt;
    shapes{nShape,1} = [G.a G.b G.c];

    
%     Ndt_SliceSpoiler = round(SliceSpoilerDur/dt);
%     Ndt_Slice = round((DwGradTsep-2*SliceSpoilerDur-DwDelay1-DwDelay2)/dt); Ndt_Slice = 2*ceil(Ndt_Slice/2);
%     Ndt_G0 = round((Ndt_GradTsep - 2*Ndt_SliceSpoiler - Ndt_Slice)/2);

%     ind_slice = logical([zeros(Ndt,1); ones(Ndt_GradTsep,1); zeros(Ndt,1)]);
%     Ndt = numel(ind_slice);
%     ind_pre180 = logical([ones(Ndt/2,1); zeros(Ndt/2,1)]); 
%     ind_pre180 = repmat(ind_pre180,[1, td1pershape]);

%     G.slice = [zeros(Ndt_G0,1); Gslicespoil*ones(Ndt_SliceSpoiler,1); Gslice*ones(Ndt_Slice,1); ...
%         Gslicespoil*ones(Ndt_SliceSpoiler,1); zeros(Ndt_G0,1)];        
%     G.slice = [zeros(Ndt_GradTsep,1)];        

%     G.a = [G.a; zeros(Ndt_GradTsep,1); G.a];
%     G.b = [G.b; zeros(Ndt_GradTsep,1); G.b];
%     G.c = [G.c; zeros(Ndt_GradTsep,1); G.c];
    G.a = [zeros(numel(G.exc_x),1); G.a; zeros(numel(G.ref_x),1); G.a; zeros(numel(G.acq_x),1)];
    G.b = [zeros(numel(G.exc_x),1); G.b; zeros(numel(G.ref_x),1); G.b; zeros(numel(G.acq_x),1)];
    G.c = [zeros(numel(G.exc_x),1); G.c; zeros(numel(G.ref_x),1); G.c; zeros(numel(G.acq_x),1)];

    if Ndt > 1
%         G.x = repmat(G.a,[1, td1pershape]).*repmat(G.xa(ind)',[Ndt, 1]) + ...
%             repmat(G.b,[1, td1pershape]).*repmat(G.xb(ind)',[Ndt, 1]) + ...
%             repmat(G.c,[1, td1pershape]).*repmat(G.xc(ind)',[Ndt, 1]);
%         G.y = repmat(G.a,[1, td1pershape]).*repmat(G.ya(ind)',[Ndt, 1]) + ...
%             repmat(G.b,[1, td1pershape]).*repmat(G.yb(ind)',[Ndt, 1]) + ...
%             repmat(G.c,[1, td1pershape]).*repmat(G.yc(ind)',[Ndt, 1]);
%         G.z = repmat(G.a,[1, td1pershape]).*repmat(G.za(ind)',[Ndt, 1]) + ...
%             repmat(G.b,[1, td1pershape]).*repmat(G.zb(ind)',[Ndt, 1]) + ...
%             repmat(G.c,[1, td1pershape]).*repmat(G.zc(ind)',[Ndt, 1]);

        G.x = repmat(G.a,[1, td1pershape]).*repmat(G.xa(ind)',[Ndt_tot, 1]) + ...
            repmat(G.b,[1, td1pershape]).*repmat(G.xb(ind)',[Ndt_tot, 1]) + ...
            repmat(G.c,[1, td1pershape]).*repmat(G.xc(ind)',[Ndt_tot, 1]);
        G.y = repmat(G.a,[1, td1pershape]).*repmat(G.ya(ind)',[Ndt_tot, 1]) + ...
            repmat(G.b,[1, td1pershape]).*repmat(G.yb(ind)',[Ndt_tot, 1]) + ...
            repmat(G.c,[1, td1pershape]).*repmat(G.yc(ind)',[Ndt_tot, 1]);
        G.z = repmat(G.a,[1, td1pershape]).*repmat(G.za(ind)',[Ndt_tot, 1]) + ...
            repmat(G.b,[1, td1pershape]).*repmat(G.zb(ind)',[Ndt_tot, 1]) + ...
            repmat(G.c,[1, td1pershape]).*repmat(G.zc(ind)',[Ndt_tot, 1]);

%         figure(2), clf, plot(dt*(1:Ndt_tot)',G.x,'r-',dt*(1:Ndt_tot)',G.y,'g-',dt*(1:Ndt_tot)',G.z,'b-'), %pause%(.1) 

%         G.x(ind_slice,:) = slicevector(1)*repmat(G.slice,[1, td1pershape]);
%         G.y(ind_slice,:) = slicevector(2)*repmat(G.slice,[1, td1pershape]);
%         G.z(ind_slice,:) = slicevector(3)*repmat(G.slice,[1, td1pershape]);
        G.x(ind_im,:) = repmat(G.im_x(ind_im),[1, td1pershape]);
        G.y(ind_im,:) = repmat(G.im_y(ind_im),[1, td1pershape]);
        G.z(ind_im,:) = repmat(G.im_z(ind_im),[1, td1pershape]);

%         G.x(ind_pre180) = -G.x(ind_pre180);
%         G.y(ind_pre180) = -G.y(ind_pre180);
%         G.z(ind_pre180) = -G.z(ind_pre180);
        G.x(ind_preref,:) = -G.x(ind_preref,:);
        G.y(ind_preref,:) = -G.y(ind_preref,:);
        G.z(ind_preref,:) = -G.z(ind_preref,:);

%         figure(1), clf, plot(dt*(1:Ndt)',G.x,'r-',dt*(1:Ndt)',G.y,'g-',dt*(1:Ndt)',G.z,'b-'), pause%(.1) 
       
        % Dephasing vector F.x in SI
        F.x = gamma*cumsum(G.x*dt,1);
        F.y = gamma*cumsum(G.y*dt,1);
        F.z = gamma*cumsum(G.z*dt,1);
        F.r = sqrt(F.x.^2 + F.y.^2 + F.z.^2); % Magnitude
%         figure(3), clf, plot(dt*(1:Ndt_tot)',F.x,'r-',dt*(1:Ndt_tot)',F.y,'g-',dt*(1:Ndt_tot)',F.z,'b-'), pause(.1) 
        plot(axh_grad,dt*(1:Ndt_tot)',G.x,'r-',dt*(1:Ndt_tot)',G.y,'g-',dt*(1:Ndt_tot)',G.z,'b-') 
        plot(axh_deph,dt*(1:Ndt_tot)',F.x/2/pi/1e6,'r-',dt*(1:Ndt_tot)',F.y/2/pi/1e6,'g-',dt*(1:Ndt_tot)',F.z/2/pi/1e6,'b-') 

        % Diffusion weighting matrix b
        bt.xx = sum(F.x.*F.x*dt,1)';
        bt.xy = sum(F.x.*F.y*dt,1)';
        bt.xz = sum(F.x.*F.z*dt,1)';
        bt.yy = sum(F.y.*F.y*dt,1)';
        bt.yz = sum(F.y.*F.z*dt,1)';
        bt.zz = sum(F.z.*F.z*dt,1)';

        b = (bt.xx + bt.yy + bt.zz); % trace

        % Save as xps
        xps.b(ind,1) = b;
        xps.bt(ind,:) = [bt.xx bt.yy bt.zz sqrt(2)*[bt.xy bt.xz bt.yz]];

        xps.gwf_x(ind,:) = G.x';
        xps.gwf_y(ind,:) = G.y';
        xps.gwf_z(ind,:) = G.z';
    end
end

% xps_old = xps;
% xps = mdm_xps_calc_btpars(xps_old);

xps_temp = xps;

% xps_temp.gwf_x = G.x'; % Maybe save full gradient modulations in the future
% xps_temp.gwf_y = G.y';
% xps_temp.gwf_z = G.z';
xps_temp.gwf_t = repmat(dt*(1:1:Ndt),[td1 1]);
%figure(1), clf, plot(xps_temp.gwf_t,xps_temp.gwf_x,'r-')

xps_temp = mdm_xps_calc_btpars(xps_temp);
xps_temp = mdm_xps_add_btomega(xps_temp,rps);

clear xps

xps.n = xps_temp.n;
xps.b = xps_temp.b;
xps.bt = xps_temp.bt;
xps.b_delta = xps_temp.b_delta;
xps.b_eta = xps_temp.b_eta;
xps.theta = xps_temp.theta;
xps.phi = xps_temp.phi;
xps.u = xps_temp.u;
xps.btomega = single(xps_temp.btomega); % Single to save disk space
xps.domega = xps_temp.domega;
xps.rmsomega = xps_temp.rmsomega;
xps.momega = xps_temp.momega;

xps.te = PVM_EchoTime*ones(td1,1);
xps.tr = PVM_RepetitionTime*ones(td1,1);

% figure(1), clf, plot((0:(size(xps.btomega,2)-1))',xps.btomega','-'), pause(1)
plot(axh_spec,linspace(0,6*rps.maxomega/2/pi,size(xps.btomega,2))',xps.btomega'/xps.domega(1),'-'), pause(1)
set([axh_grad; axh_deph; axh_spec],'Box','off','TickDir','out','LineWidth',1,'FontSize',8)
set([axh_grad; axh_deph],'XLim',dt*Ndt_tot*[-.05 1.05])
ylabel(axh_grad,'g / Tm^{-1}')
ylabel(axh_deph,'q / 10^6 m^{-1}')
ylabel(axh_spec,'b(\omega) / s^2m^{-2}')
xlabel(axh_deph,'t / s')
xlabel(axh_spec,'\omega / 2\pi Hz (Voigt-like xx, yy, zz, xy, xz, yz)')
set(axh_grad,'XTickLabel',[])
set(axh_spec,'XLim',6*rps.maxomega/2/pi*[-.05 1.05],'XTick',[0 rps.maxomega/2/pi])
papersize = 2*[8.3 5.9];
set(fh_shapes, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
eval(['print ' fileparts(xps_fn) '/waveforms -loose -dpdf'])


% xps.te = EchoTime*ones(xps.n,1);

% if (~isempty(DwNDiffExp) && strcmp(MultipleEchos,'yes'))
% 
%     Nte = numel(EffectiveEchoTime);
%     xps_cell = cell(Nte,1);
%     for nte = 1:Nte
%         xps_cell{nte} = xps;
%         xps_cell{nte}.te =EffectiveEchoTime(nte)*ones(xps.n,1);
%     end
%     xps = mdm_xps_merge(xps_cell);
%     
% 
%     EpiAcquistionTime = ReadPVParam(data_path_pv, 'EpiAcquistionTime')*1e-3 ;
%     PVM_Matrix = ReadPVParam(data_path_pv, 'PVM_Matrix');
% 
%     %xps.dte = repmat(EpiAcquistionTime*fliplr(-1/2:1/(PVM_Matrix(1,2)-1):1/2),[xps.n 1]);
%     xps.dte = repmat(EpiAcquistionTime*(-1/2:1/(PVM_Matrix(1,2)-1):1/2),[xps.n 1]);
% 
% end

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

if (~isempty(xps_fn))
    save(xps_fn,'xps');
%     save(fullfile(fileparts(xps_fn),'shapes'),'shapes'); %Uncomment to save gradient shapes to file
end


