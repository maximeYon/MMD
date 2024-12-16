function xps = mdm_bruker_my_DOR_R12_RARE_PV360_acqus2xps(data_path, xps_fn, rps)
% function xps = mdm_bruker_my_ll_R1_R2_DTD_MSME_PV360_acqus2xps(data_path, xps_fn)
%
% Calculation of b-tensors for gradient waveforms
% from Topgaard, Phys. Chem. Chem. Phys., 2016.
% http://dx.doi.org/10.1039/c5cp07251d
% Written by Maxime Yon (2021)
%
% data_path: directory for gradient text files
% xps_fn (optional) filename of xps

if (nargin < 2), xps_fn = []; end


%% Get parameters
data_path_pv = [data_path filesep];
DwNShapes = ReadPV360Param(data_path_pv, 'DwNShapes') ;
DiffusionExperimentList = ReadPV360Param(data_path_pv, 'DiffusionExperimentList') ;
ShapeNumberList = ReadPV360Param(data_path_pv, 'ShapeNumberList') ;
EffectiveTEList = ReadPV360Param(data_path_pv, 'EffectiveTEList')*1e-3 ;
EffectiveTRList = (ReadPV360Param(data_path_pv, 'EffectiveTRList'))'*1e-3 ;
PVM_GradCalConst = ReadPV360Param(data_path_pv, 'PVM_GradCalConst') ;
ExcPulse1=ReadPV360Param(data_path_pv, 'ExcPulse1') ;
ExcPulse1_Rpfac =  split(ExcPulse1,',');
ExcPulse1_Rpfac = str2double(ExcPulse1_Rpfac{10})/100;
gradAmpl=ReadPV360Param(data_path_pv, 'ACQ_gradient_amplitude') ;
Delays=ReadPV360Param(data_path_pv, 'D') ; % in seconde
Pulses=ReadPV360Param(data_path_pv, 'P') ; % in useconde
RefoCrusher=ReadPV360Param(data_path_pv, 'RefoCrusherYesNo') ;
%% RAREst parameters
%param
PhaseEncList=ReadPV360Param(data_path_pv, 'ACQ_spatial_phase_1');
NSpinEcho=ReadPV360Param(data_path_pv, 'PVM_EncCentralStep1');
isNUS = ReadPV360Param(data_path_pv, 'CSYesNo');
if strcmp(isNUS,'yes')
   NSpinEcho = 1; 
end
% NSpinEcho = find(PhaseEncList==0);
Method_name = ReadPV360Param(data_path_pv, 'Method');
VFA=ReadPV360Param(data_path_pv, 'VFAYesNo') ;
if strcmp(VFA,'yes')==1
    VFA_FAList=ReadPV360Param(data_path_pv, 'VFA_FAList') ;
    VFA_PESkips=ReadPV360Param(data_path_pv, 'VFA_PESkips') ;
    if strcmp(Method_name,'<user:my_dor12_rare_final>')==1 
    RareFactor=ReadPV360Param(data_path_pv, 'MY_RareFactor') ;    
    else
    RareFactor=ReadPV360Param(data_path_pv, 'PVM_RareFactor') ;
    end
    NSpinEcho = NSpinEcho + VFA_PESkips;
end

BIRRefoc=ReadPV360Param(data_path_pv, 'BIRRefocYesNo') ;
if strcmp(BIRRefoc,'yes')==1
    EchoPulse1=ReadPV360Param(data_path_pv, 'EchoPulse1') ;
    EchoPulse2=ReadPV360Param(data_path_pv, 'EchoPulse2') ;
end
%% Retrieve baselevel parameters
% grad val
g0 = gradAmpl(0+1);
g2 = gradAmpl(2+1);
g3 = gradAmpl(3+1);
g4 = gradAmpl(4+1);
g5 = gradAmpl(5+1);
g6 = 0; %force center of Kspace
g7 = gradAmpl(7+1);
g8 = gradAmpl(8+1);
g12 = gradAmpl(12+1);
if strcmp(RefoCrusher,'yes')==1
    g22 = gradAmpl(22+1);
end
% g0 =0;g2 =0;g3 =0;g4 =0;g5 =0;g7 =0;g8 =0;g12 =0;g22=0;
% pulse duration
p0 = Pulses(0+1)*1e-6;
p1 = Pulses(1+1)*1e-6;
if strcmp(BIRRefoc,'yes')==1
    p3 = Pulses(3+1)*1e-6;
    p4 = Pulses(4+1)*1e-6;
end
%delays
d1 = Delays(1+1);
d2 = Delays(2+1);
d3 = Delays(3+1);
d5 = Delays(5+1);
d6 = Delays(6+1);
d11 = Delays(11+1);
d12 = Delays(12+1);
denab = d3;
dAq = ReadPV360Param(data_path_pv, 'PVM_AcquisitionTime')*1e-3;
if strcmp(RefoCrusher,'yes')==1
    d22 = Delays(22+1);
end

AdditionalTESec=ReadPV360Param(data_path_pv, 'AdditionalTESec') ; % in seconde, on each side of echo
DwD1 = ReadPV360Param(data_path_pv, 'DwD1');
DwD4 = ReadPV360Param(data_path_pv, 'DwD4');

%% Check Repetition time
d0 = Delays(0+1);
d9 = Delays(9+1);
d10 = Delays(10+1);
AddTRList = (ReadPV360Param(data_path_pv, 'AdditionalTRSec'));
u10 = 0.000010;
EffecectiveTRcalc = d0+d10+AddTRList+u10+d9+d3+(p0*(1-ExcPulse1_Rpfac));

% figure()
% hold on
% plot(EffectiveTRList,'r');
% plot(EffecectiveTRcalc,'b');

%% Compute shapes
dt = 0.5*8e-6; % 4 useconde
% dt = 0.4*10e-6; % 4 useconde
Ndt_DwD4 = round(DwD4/dt);

for nShape = 1:DwNShapes
    ShapeStru.(['DWGS' num2str(nShape-1) 'R']) = ReadPV360Param(data_path_pv,['DWGS' num2str(nShape-1) 'R']);
    ShapeStru.(['DWGS' num2str(nShape-1) 'P']) = ReadPV360Param(data_path_pv,['DWGS' num2str(nShape-1) 'P']);
    ShapeStru.(['DWGS' num2str(nShape-1) 'S']) = ReadPV360Param(data_path_pv,['DWGS' num2str(nShape-1) 'S']);
    
    % Avoid error if shape is one points
    if numel(ShapeStru.(['DWGS' num2str(nShape-1) 'R'])) == 1
        Ndt = round(DwD4/dt);
        ShapeStru.(['DWGS' num2str(nShape-1) 'R']) = ShapeStru.(['DWGS' num2str(nShape-1) 'R'])*ones(Ndt,1);
        ShapeStru.(['DWGS' num2str(nShape-1) 'P']) = ShapeStru.(['DWGS' num2str(nShape-1) 'P'])*ones(Ndt,1);
        ShapeStru.(['DWGS' num2str(nShape-1) 'S']) = ShapeStru.(['DWGS' num2str(nShape-1) 'S'])*ones(Ndt,1);
    end
end
Size_shape = size(ShapeStru.DWGS0R,2);
%% Read gradient values
% G.xa etc maps gradient modulations (a,b,c) to channels (x,y,z)
% Normalized from -1 to +1
G.xa = ReadPV360Param(data_path_pv, 'DwGxa')' ;
G.xb = ReadPV360Param(data_path_pv, 'DwGxb')' ;
G.xc = ReadPV360Param(data_path_pv, 'DwGxc')' ;
G.ya = ReadPV360Param(data_path_pv, 'DwGya')' ;
G.yb = ReadPV360Param(data_path_pv, 'DwGyb')' ;
G.yc = ReadPV360Param(data_path_pv, 'DwGyc')' ;
G.za = ReadPV360Param(data_path_pv, 'DwGza')' ;
G.zb = ReadPV360Param(data_path_pv, 'DwGzb')' ;
G.zc = ReadPV360Param(data_path_pv, 'DwGzc')' ;

td1 = numel(G.xa);

% Convert ramps r.ax to uvec axis of q-vector cone
gpas.avecnorm = sqrt(G.xa.^2+G.ya.^2+G.za.^2);
gpas.bvecnorm = sqrt(G.xb.^2+G.yb.^2+G.zb.^2);
gpas.cvecnorm = sqrt(G.xc.^2+G.yc.^2+G.zc.^2);
gpas.avec = [G.xa G.ya G.za]./repmat(gpas.avecnorm,[1 3]);
gpas.bvec = [G.xb G.yb G.zb]./repmat(gpas.bvecnorm,[1 3]);
gpas.cvec = [G.xc G.yc G.zc]./repmat(gpas.cvecnorm,[1 3]);
gpas.abcrossvec = msf_crossprod_nx3vectors(gpas.avec,gpas.bvec);
uvec = zeros(td1,3);
ind = find(gpas.avecnorm>0);
uvec(ind,:) = gpas.abcrossvec(ind,:);
ind = find(gpas.cvecnorm>0);
uvec(ind,:) = gpas.cvec(ind,:);

gamma = 26.75e7;
field_names_G = fieldnames(G);
for ind_field = 1:size(field_names_G,1)
    G.(field_names_G{ind_field}) = G.(field_names_G{ind_field})/100*PVM_GradCalConst*2*pi/gamma*1e3;
end

xps.n = td1;
xps.b = zeros(td1,1);
xps.bt = zeros(td1,6);
xps.u = uvec;
xps.theta = acos(uvec(:,3));
xps.phi = atan2(uvec(:,2),uvec(:,1));

%% Compute delay with dt
Ndt_p0half = round(p0*ExcPulse1_Rpfac/dt);
Ndt_p1half = round(p1/2/dt);
if strcmp(BIRRefoc,'yes')==1
 Ndt_p3 = round(p3/dt);  
 Ndt_p4 = round(p4/dt); 
end

Ndt_d1 = round(d1/dt);
Ndt_d2 = round(d2/dt);
Ndt_d3 = round(d3/dt);
Ndt_d5 = round(d5/dt);
Ndt_d6 = round(d6/dt);
Ndt_d11 = round(d11/dt);
Ndt_d12 = round(d12/dt);
Ndt_Aq = round(dAq/dt);
Ndt_AqDiv2 = round(Ndt_Aq/2);
Ndt_denab = round(denab/dt);

Ndt_AddTE = round(AdditionalTESec/dt);
Ndt_DwD1 = round(DwD1/dt);
if strcmp(RefoCrusher,'yes')==1
    Ndt_d22 = round(d22/dt);
end

%% Compute main sequence part timing
if strcmp(BIRRefoc,'yes')==1
    Ndt_preref = Ndt_p0half + Ndt_d3 + Ndt_d12 + Ndt_d3 +...
    Ndt_AddTE + (Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1).*DiffusionExperimentList  +...
    Ndt_d3 + Ndt_p3 + Ndt_p4;
else
    Ndt_preref = Ndt_p0half + Ndt_d3 + Ndt_d12 + Ndt_d3 +...
    Ndt_AddTE + (Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1).*DiffusionExperimentList  +...
    Ndt_d3 + Ndt_p1half;
end
if strcmp(RefoCrusher,'yes')==1
    Ndt_preref = Ndt_preref + Ndt_d22;
end

if strcmp(BIRRefoc,'yes')==1
    Ndt_postref = Ndt_p2 + Ndt_p3 + Ndt_p4 + Ndt_d3 +...
    (Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1).*DiffusionExperimentList + Ndt_AddTE +  ...
    Ndt_d3 + Ndt_d12 + Ndt_d3 + Ndt_d11 + Ndt_d2 + Ndt_d1 +...
    (Ndt_d3 + 2*Ndt_p1half + Ndt_d5 + Ndt_d6 + Ndt_denab + Ndt_Aq + Ndt_d6 + Ndt_d5).*(NSpinEcho-1) +... % loop
    Ndt_d3 + 2*Ndt_p1half  + Ndt_d5 + Ndt_d6 + Ndt_denab+ Ndt_AqDiv2;
else
Ndt_postref = Ndt_p1half + Ndt_d3 +...
    (Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1).*DiffusionExperimentList + Ndt_AddTE +  ...
    Ndt_d3 + Ndt_d12 + Ndt_d3 + Ndt_d11 + Ndt_d2 + Ndt_d1 +...
    (Ndt_d3 + 2*Ndt_p1half + Ndt_d5 + Ndt_d6 + Ndt_denab + Ndt_Aq + Ndt_d6 + Ndt_d5).*(NSpinEcho-1) +... % loop
    Ndt_d3 + 2*Ndt_p1half  + Ndt_d5 + Ndt_d6 + Ndt_denab+ Ndt_AqDiv2;
end
if strcmp(RefoCrusher,'yes')==1
    Ndt_postref = Ndt_postref + Ndt_d22;
end
Ndt_tot = Ndt_preref + Ndt_postref;

%% Create sequences part for b-calculation
Seq_part(1,:) = Ndt_p0half + Ndt_d3 + Ndt_d12 + Ndt_d3 +Ndt_AddTE;
Seq_part(2,:) = Seq_part(1,:) + (Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1).*DiffusionExperimentList;
if strcmp(BIRRefoc,'yes')==1
Seq_part(3,:) = Seq_part(2,:) + Ndt_d3 + Ndt_p3 + Ndt_p4 + Ndt_p3 + Ndt_p4 + Ndt_d3;    
else
Seq_part(3,:) = Seq_part(2,:) + Ndt_d3 + Ndt_p1half + Ndt_p1half + Ndt_d3;
end
if strcmp(RefoCrusher,'yes')==1
    Seq_part(3,:) = Seq_part(3,:) + 2* Ndt_d22;
end
Seq_part(4,:) = Seq_part(3,:) + (Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1).*DiffusionExperimentList + Ndt_AddTE;
Seq_part(5,:) = Seq_part(4,:) + Ndt_d3 + Ndt_d12 + Ndt_d3 + Ndt_d11 + Ndt_d2 + Ndt_d1 +...
    (Ndt_d3 + 2*Ndt_p1half + Ndt_d5 + Ndt_d6 + Ndt_denab + Ndt_Aq + Ndt_d6 + Ndt_d5).*(NSpinEcho-1) +... % loop
    Ndt_d3 + 2*Ndt_p1half  + Ndt_d5 + Ndt_d6 + Ndt_denab+ Ndt_AqDiv2;

xps.gwf_x = zeros(td1,max(Ndt_tot)); %Use max time step
xps.gwf_y = zeros(td1,max(Ndt_tot));
xps.gwf_z = zeros(td1,max(Ndt_tot));

%% Special TE calculation for VFA
if strcmp(VFA,'yes')==1
Spin_echo_dur_right = Ndt_preref*dt;
Spin_echo_dur_left = (Ndt_postref - (Ndt_d2 + Ndt_d1 + (Ndt_d3 + 2*Ndt_p1half + Ndt_d5 + Ndt_d6 + Ndt_denab + Ndt_Aq + Ndt_d6 + Ndt_d5).*(NSpinEcho-1) +... % loop
    Ndt_d3 + 2*Ndt_p1half  + Ndt_d5 + Ndt_d6 + Ndt_denab+ Ndt_AqDiv2))*dt;
Spin_echo_dur = Spin_echo_dur_right + Spin_echo_dur_left; % in s
Inter_RARE_echo_delay = (Ndt_d3 + 2*Ndt_p1half + Ndt_d5 + Ndt_d6 + Ndt_denab + Ndt_Aq + Ndt_d6 + Ndt_d5)*dt; % in s
VFA_FAList = VFA_FAList./180.*pi();
Number_of_echoes = VFA_PESkips + RareFactor;
VFA_FAList = VFA_FAList(1,1:Number_of_echoes);
for ind_TE = 1:size(Spin_echo_dur,2)
    VFA_TE(ind_TE,:) = my_calc_TE_vfa(Number_of_echoes,Spin_echo_dur(ind_TE),Inter_RARE_echo_delay,VFA_FAList);
end
VFA_TE = VFA_TE(:,NSpinEcho);
end

%% Check echo time calculations
NdtTE = Ndt_tot*dt;
% figure()
% hold on
% plot(NdtTE,'b')
% plot(EffectiveTEList,'r')
% plot(VFA_TE,'k')

unit2si = 1/100*PVM_GradCalConst*2*pi/gamma*1e3;
% Recompute g0
g0_exc = - g12*Ndt_d12/(Ndt_p0half+Ndt_d3/2);
%% Asociate delays with gradients
G.exc_x = [0*ones(Ndt_p0half,1); 0*ones(Ndt_d3,1); 0*ones(Ndt_d12,1); 0*ones(Ndt_d3,1)]*unit2si;
G.exc_y = [0*ones(Ndt_p0half,1); 0*ones(Ndt_d3,1); 0*ones(Ndt_d12,1); 0*ones(Ndt_d3,1)]*unit2si;
G.exc_z = [g0_exc*ones(Ndt_p0half,1); g0_exc/2*ones(Ndt_d3,1); g12*ones(Ndt_d12,1); 0*ones(Ndt_d3,1)]*unit2si;

if strcmp(BIRRefoc,'yes')==1
    if strcmp(RefoCrusher,'yes')==1
    G.ref_x = [0*ones(Ndt_d22,1);0*ones(Ndt_d3,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1);0*ones(Ndt_d22,1); 0*ones(Ndt_d3,1);]*unit2si;
    G.ref_y = [0*ones(Ndt_d22,1);0*ones(Ndt_d3,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1);0*ones(Ndt_d22,1); 0*ones(Ndt_d3,1);]*unit2si;
    G.ref_z = [g22*ones(Ndt_d22,1);0*ones(Ndt_d3,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1);g22*ones(Ndt_d22,1); 0*ones(Ndt_d3,1);]*unit2si;
else
    G.ref_x = [0*ones(Ndt_d3,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1); 0*ones(Ndt_d3,1);]*unit2si;
    G.ref_y = [0*ones(Ndt_d3,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1); 0*ones(Ndt_d3,1);]*unit2si;
    G.ref_z = [0*ones(Ndt_d3,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1);0*ones(Ndt_p3,1); 0*ones(Ndt_p4,1); 0*ones(Ndt_d3,1);]*unit2si;
end
else
if strcmp(RefoCrusher,'yes')==1
    G.ref_x = [0*ones(Ndt_d22,1);0*ones(Ndt_d3,1);0*ones(Ndt_p1half,1); 0*ones(Ndt_p1half,1);0*ones(Ndt_d22,1); 0*ones(Ndt_d3,1);]*unit2si;
    G.ref_y = [0*ones(Ndt_d22,1);0*ones(Ndt_d3,1);0*ones(Ndt_p1half,1); 0*ones(Ndt_p1half,1);0*ones(Ndt_d22,1); 0*ones(Ndt_d3,1);]*unit2si;
    G.ref_z = [g22*ones(Ndt_d22,1);g0/2*ones(Ndt_d3,1);g0*ones(Ndt_p1half,1); g0*ones(Ndt_p1half,1);g22*ones(Ndt_d22,1); g0/2*ones(Ndt_d3,1);]*unit2si;
else
    G.ref_x = [0*ones(Ndt_d3,1);0*ones(Ndt_p1half,1); 0*ones(Ndt_p1half,1); 0*ones(Ndt_d3,1);]*unit2si;
    G.ref_y = [0*ones(Ndt_d3,1);0*ones(Ndt_p1half,1); 0*ones(Ndt_p1half,1); 0*ones(Ndt_d3,1);]*unit2si;
    G.ref_z = [g0/2*ones(Ndt_d3,1);g0*ones(Ndt_p1half,1); g0*ones(Ndt_p1half,1); g0/2*ones(Ndt_d3,1);]*unit2si;
end
end

G.acq_x = [0*ones(Ndt_d3,1);0*ones(Ndt_d12,1);0*ones(Ndt_d3,1); 0*ones(Ndt_d11,1); g2*ones(Ndt_d2,1);0*ones(Ndt_d1,1)]*unit2si; % decrease grad ramp during d1;
G.acq_y = [0*ones(Ndt_d3,1);0*ones(Ndt_d12,1);0*ones(Ndt_d3,1); 0*ones(Ndt_d11,1); 0*ones(Ndt_d2,1);0*ones(Ndt_d1,1)]*unit2si;
G.acq_z = [0*ones(Ndt_d3,1);0*ones(Ndt_d12,1);0*ones(Ndt_d3,1); 0*ones(Ndt_d11,1); g7*ones(Ndt_d2,1);0*ones(Ndt_d1,1)]*unit2si;

for ind_echo = 1:NSpinEcho-1
    if  mod(ind_echo,2) % Odd
        SignBpulse = 1; SignApulse = -1;
    else
        SignBpulse = -1; SignApulse = 1;
    end
    
    if strcmp(VFA,'yes')==1
        if ind_echo<VFA_PESkips+1
            g4_eff = 0;
        else
            g4_eff = g4*PhaseEncList(1,ind_echo-VFA_PESkips);
        end
    else
        g4_eff = g4*PhaseEncList(1,ind_echo);
    end
    
    G.acq_x = cat(1,G.acq_x,[SignBpulse*0*ones(Ndt_d3,1);SignBpulse*0*ones(Ndt_p1half,1);SignApulse*0*ones(Ndt_p1half,1);SignApulse*g8*ones(Ndt_d5,1);SignApulse*0*ones(Ndt_d6,1);SignApulse*g5*ones(Ndt_denab,1);SignApulse*g5*ones(Ndt_Aq,1);SignApulse*0*ones(Ndt_d6,1);SignApulse*g8*ones(Ndt_d5,1)]*unit2si);
    G.acq_y = cat(1,G.acq_y,[SignBpulse*0*ones(Ndt_d3,1);SignBpulse*0*ones(Ndt_p1half,1);SignApulse*0*ones(Ndt_p1half,1);SignApulse*g4_eff*ones(Ndt_d5,1);SignApulse*0*ones(Ndt_d6,1);SignApulse*0*ones(Ndt_denab,1);SignApulse*0*ones(Ndt_Aq,1);SignApulse*0*ones(Ndt_d6,1);SignApulse*g4_eff*ones(Ndt_d5,1)]*unit2si);
    G.acq_z = cat(1,G.acq_z,[SignBpulse*0*ones(Ndt_d3,1);SignBpulse*g0*ones(Ndt_p1half,1);SignApulse*g0*ones(Ndt_p1half,1);SignApulse*g3-g6*ones(Ndt_d5,1);SignApulse*0*ones(Ndt_d6,1);SignApulse*0*ones(Ndt_denab,1);SignApulse*0*ones(Ndt_Aq,1);SignApulse*0*ones(Ndt_d6,1);SignApulse*g3+g6*ones(Ndt_d5,1)]*unit2si);
end
ind_echo = ind_echo+1;
if isempty(ind_echo); ind_echo =1; end

if  mod(ind_echo,2) % Odd
    SignBpulse = 1; SignApulse = -1;
else
    SignBpulse = -1; SignApulse = 1;
end

if strcmp(VFA,'yes')==1
    if ind_echo<VFA_PESkips+1
        g4_eff = 0;
    else
        g4_eff = g4*PhaseEncList(1,ind_echo-VFA_PESkips);
    end
else
    g4_eff = g4*PhaseEncList(1,ind_echo);
end

G.acq_x = cat(1,G.acq_x,[SignBpulse*0*ones(Ndt_d3,1);SignBpulse*0*ones(Ndt_p1half,1);SignApulse*0*ones(Ndt_p1half,1);SignApulse*g8*ones(Ndt_d5,1);SignApulse*0*ones(Ndt_d6,1);SignApulse*g5*ones(Ndt_denab,1);SignApulse*g5*ones(Ndt_AqDiv2,1)]*unit2si);
G.acq_y = cat(1,G.acq_y,[SignBpulse*0*ones(Ndt_d3,1);SignBpulse*0*ones(Ndt_p1half,1);SignApulse*0*ones(Ndt_p1half,1);SignApulse*g4_eff*ones(Ndt_d5,1);SignApulse*0*ones(Ndt_d6,1);SignApulse*0*ones(Ndt_denab,1);SignApulse*0*ones(Ndt_AqDiv2,1)]*unit2si);
G.acq_z = cat(1,G.acq_z,[SignBpulse*0*ones(Ndt_d3,1);SignBpulse*g0*ones(Ndt_p1half,1);SignApulse*g0*ones(Ndt_p1half,1);SignApulse*g3-g6*ones(Ndt_d5,1);SignApulse*0*ones(Ndt_d6,1);SignApulse*0*ones(Ndt_denab,1);SignApulse*0*ones(Ndt_AqDiv2,1)]*unit2si);

%% combine imaging gradients
Ndt_DwGradDur = Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1;
ShapeNdtPos = zeros(td1,2);
for ind_image = 1:td1
    if DiffusionExperimentList(ind_image)==1
        ShapeNdtPos(ind_image,:) = [size(G.exc_x,1)+Ndt_AddTE(ind_image)+Ndt_DwD1+1 size(G.exc_x,1)+Ndt_AddTE(ind_image)+Ndt_DwGradDur(ind_image)+size(G.ref_x,1)+Ndt_DwD1+1];
        G.im_x{ind_image} = [G.exc_x; zeros(Ndt_AddTE(ind_image),1); zeros(Ndt_DwGradDur(ind_image),1); G.ref_x; zeros(Ndt_DwGradDur(ind_image),1); zeros(Ndt_AddTE(ind_image),1); G.acq_x];
        G.im_y{ind_image} = [G.exc_y; zeros(Ndt_AddTE(ind_image),1); zeros(Ndt_DwGradDur(ind_image),1); G.ref_y; zeros(Ndt_DwGradDur(ind_image),1); zeros(Ndt_AddTE(ind_image),1); G.acq_y];
        G.im_z{ind_image} = [G.exc_z; zeros(Ndt_AddTE(ind_image),1); zeros(Ndt_DwGradDur(ind_image),1); G.ref_z; zeros(Ndt_DwGradDur(ind_image),1); zeros(Ndt_AddTE(ind_image),1); G.acq_z];
    else
        G.im_x{ind_image} = [G.exc_x; zeros(Ndt_AddTE(ind_image),1); G.ref_x; zeros(Ndt_AddTE(ind_image),1); G.acq_x];
        G.im_y{ind_image} = [G.exc_y; zeros(Ndt_AddTE(ind_image),1); G.ref_y; zeros(Ndt_AddTE(ind_image),1); G.acq_y];
        G.im_z{ind_image} = [G.exc_z; zeros(Ndt_AddTE(ind_image),1); G.ref_z; zeros(Ndt_AddTE(ind_image),1); G.acq_z];
    end
end

fh_shapes = figure(1); clf
axh_grad = axes(fh_shapes,'position',[.1 .73 .85 .25]);
axh_deph = axes(fh_shapes,'position',[.1 .45 .85 .25]);
axh_spec = axes(fh_shapes,'position',[.1 .1 .85 .25]);
hold(axh_grad,'on')
hold(axh_deph,'on')
hold(axh_spec,'on')

%% Add shapes and orientation
for ind_image = 1:td1
    nShape = ShapeNumberList(ind_image);
    if DiffusionExperimentList(ind_image)==1
        % retrieve current shape
        G.a = ShapeStru.(['DWGS' num2str(nShape) 'R']);
        G.b = ShapeStru.(['DWGS' num2str(nShape) 'P']);
        G.c = ShapeStru.(['DWGS' num2str(nShape) 'S']);
        
        % Interp shape
        G.a = interp1(linspace(0,1,Size_shape)',ShapeStru.(['DWGS' num2str(nShape) 'R']),linspace(0,1,Ndt_DwD4(ind_image))','makima');
        G.b = interp1(linspace(0,1,Size_shape)',ShapeStru.(['DWGS' num2str(nShape) 'P']),linspace(0,1,Ndt_DwD4(ind_image))','makima');
        G.c = interp1(linspace(0,1,Size_shape)',ShapeStru.(['DWGS' num2str(nShape) 'S']),linspace(0,1,Ndt_DwD4(ind_image))','makima');
        
        % Combine shape with coef xa,xb,... zc
        GshapeX = G.a.*G.xa(ind_image) + G.b.*G.xb(ind_image) + G.c.*G.xc(ind_image);
        GshapeY = G.a.*G.ya(ind_image) + G.b.*G.yb(ind_image) + G.c.*G.yc(ind_image);
        GshapeZ = G.a.*G.za(ind_image) + G.b.*G.zb(ind_image) + G.c.*G.zc(ind_image);
        % Add shape to imaging gradients
        G.x = G.im_x{ind_image};
        G.x(ShapeNdtPos(ind_image,1):ShapeNdtPos(ind_image,1)+size(G.a,1)-1) = GshapeX;
        G.x(ShapeNdtPos(ind_image,2):ShapeNdtPos(ind_image,2)+size(G.a,1)-1) = GshapeX;
        G.y = G.im_y{ind_image};
        G.y(ShapeNdtPos(ind_image,1):ShapeNdtPos(ind_image,1)+size(G.b,1)-1) = GshapeY;
        G.y(ShapeNdtPos(ind_image,2):ShapeNdtPos(ind_image,2)+size(G.b,1)-1) = GshapeY;
        G.z = G.im_z{ind_image};
        G.z(ShapeNdtPos(ind_image,1):ShapeNdtPos(ind_image,1)+size(G.c,1)-1) = GshapeZ;
        G.z(ShapeNdtPos(ind_image,2):ShapeNdtPos(ind_image,2)+size(G.c,1)-1) = GshapeZ;
    else
        G.x = G.im_x{ind_image};
        G.y = G.im_y{ind_image};
        G.z = G.im_z{ind_image};
    end
    
    G.x(1:Ndt_preref(ind_image)) = -G.x(1:Ndt_preref(ind_image));
    G.y(1:Ndt_preref(ind_image)) = -G.y(1:Ndt_preref(ind_image));
    G.z(1:Ndt_preref(ind_image)) = -G.z(1:Ndt_preref(ind_image));
    
    %% Check cumsum
%         CumSumX = cumsum(G.x); CumSumY = cumsum(G.y); CumSumZ = cumsum(G.z);
%         figure(2)
%         clf
%         subplot(2,1,1)
%         hold on
%         plot(CumSumX,'r')
%         plot(CumSumY,'b')
%         plot(CumSumZ,'k')
%         subplot(2,1,2)
%         hold on
%         plot(G.x,'r')
%         plot(G.y,'b')
%         plot(G.z,'k')
%         legend('X','Y','Z')
    
    % Dephasing vector F.x in SI
    F.x = gamma*cumsum(G.x*dt,1);
    F.y = gamma*cumsum(G.y*dt,1);
    F.z = gamma*cumsum(G.z*dt,1);
    F.r = sqrt(F.x.^2 + F.y.^2 + F.z.^2); % Magnitude
    
    plot(axh_grad,dt*(1:Ndt_tot(1,ind_image))',G.x,'r-',dt*(1:Ndt_tot(1,ind_image))',G.y,'g-',dt*(1:Ndt_tot(1,ind_image))',G.z,'b-')
    plot(axh_deph,dt*(1:Ndt_tot(1,ind_image))',F.x/2/pi/1e6,'r-',dt*(1:Ndt_tot(1,ind_image))',F.y/2/pi/1e6,'g-',dt*(1:Ndt_tot(1,ind_image))',F.z/2/pi/1e6,'b-')
    
    % Diffusion weighting matrix b
    bt.xx = sum(F.x.*F.x*dt,1)';
    bt.xy = sum(F.x.*F.y*dt,1)';
    bt.xz = sum(F.x.*F.z*dt,1)';
    bt.yy = sum(F.y.*F.y*dt,1)';
    bt.yz = sum(F.y.*F.z*dt,1)';
    bt.zz = sum(F.z.*F.z*dt,1)';
    
    b = (bt.xx + bt.yy + bt.zz); % trace
    
    % Save as xps
    xps.b(ind_image,1) = b;
    xps.bt(ind_image,:) = [bt.xx bt.yy bt.zz sqrt(2)*[bt.xy bt.xz bt.yz]];
    
    xps.gwf_x(ind_image,1:size(G.x,1)) = G.x';
    xps.gwf_y(ind_image,1:size(G.y,1)) = G.y';
    xps.gwf_z(ind_image,1:size(G.z,1)) = G.z';
    
    %     drawnow
    %     pause(1)
    %     cla(axh_grad)
    %     cla(axh_deph)
end

xps_temp = xps;

xps_temp.gwf_t = repmat(dt*(1:1:max(Ndt_DwD4)),[td1 1]);
xps_temp = mdm_xps_calc_btpars(xps_temp);
% xps_temp = mdm_xps_add_btomega(xps_temp,rps);
xps_temp = my_mdm_xps_add_btomega(xps_temp,rps);
% xps_temp = my_mdm_xps_add_btomega_part(xps_temp,rps,Seq_part);

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

if strcmp(VFA,'yes')==1
    xps.te = VFA_TE;
else
%     xps.te = EffectiveTEList;
    xps.te = NdtTE';
end
% xps.tr = EffectiveTRList;
xps.tr = EffecectiveTRcalc';

% figure(1), clf, plot((0:(size(xps.btomega,2)-1))',xps.btomega','-'), pause(1)
plot(axh_spec,linspace(0,6*rps.maxomega/2/pi,size(xps.btomega,2))',xps.btomega'/xps.domega(1),'-'), pause(1)
set([axh_grad; axh_deph; axh_spec],'Box','off','TickDir','out','LineWidth',1,'FontSize',8)
set([axh_grad; axh_deph],'XLim',dt*max(Ndt_tot)*[-.05 1.05])
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

if (~isempty(xps_fn))
    save(xps_fn,'xps');
    %     save(fullfile(fileparts(xps_fn),'shapes'),'shapes'); %Uncomment to save gradient shapes to file
end


