function xps = mdm_bruker_my_DOR_R12_EPI_PV6_acqus2xps(data_path, xps_fn, rps)
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
ExcPulse1=ReadPV360Param(data_path_pv, 'ExcPul') ;
ExcPulse1_Rpfac =  split(ExcPulse1,',');
ExcPulse1_Rpfac = str2double(ExcPulse1_Rpfac{10})/100;
gradAmpl=ReadPV360Param(data_path_pv, 'ACQ_gradient_amplitude') ;
Delays=ReadPV360Param(data_path_pv, 'D') ; % in seconde
Pulses=ReadPV360Param(data_path_pv, 'P') ; % in useconde
%% EPI parameters
%Delay
EpiD4=ReadPV360Param(data_path_pv, 'PVM_EpiD4'); % PVM_EpiDephaseTime
EpiD5=ReadPV360Param(data_path_pv, 'PVM_EpiD5'); % PVM_EpiSwitchTime
EpiD11=ReadPV360Param(data_path_pv, 'PVM_EpiD11'); % PVM_EpiPlateau
% EpiNEchoes=ReadPV360Param(data_path_pv, 'PVM_EpiNEchoes') ;
% Gradients
list_shape = [1 3 18 14 9 11 8 2 4 13 12 7 10 5 6];
for indshape = list_shape
    EpiShape.(['s' num2str(indshape)])=ReadPV360Param(data_path_pv, ['EpiShape' num2str(indshape)]) ;%
end

%% Retrieve baselevel parameters
% grad val
g0 = gradAmpl(0+1);
g1 = gradAmpl(1+1);
g3 = gradAmpl(3+1);
% g0 =0; g1 =0; g3 =0;
% pulse duration
p0 = Pulses(0+1)*1e-6;
p1 = Pulses(1+1)*1e-6;
%delays
d1 = Delays(1+1);
d2 = Delays(2+1);
d3 = Delays(3+1);
d4 = Delays(4+1);
d7 = Delays(7+1);
d8 = Delays(8+1);
d9 = Delays(9+1);
d11 = Delays(11+1); % Eddy current resting delay

AdditionalTESec=ReadPV360Param(data_path_pv, 'AdditionalTESec') ; % in seconde, on each side of echo
DwD1 = ReadPV360Param(data_path_pv, 'DwD1');
DwD4 = ReadPV360Param(data_path_pv, 'DwD4');

%% Check Repetition time of multislice acquisition
Nslices=ReadPV360Param(data_path_pv, 'PVM_SPackArrNSlices') ;
AddTRList = (ReadPV360Param(data_path_pv, 'AdditionalTRSec'));
if Nslices >1
    SliceTime = (ReadPV360Param(data_path_pv, 'SliceTime'));
    EffectiveTEList = (ReadPV360Param(data_path_pv, 'EffectiveTEList'));
    EffectiveTEList = permute(repmat(EffectiveTEList,1,Nslices),[2 1]); EffectiveTEList = EffectiveTEList(:);
    
    for ind_rep = 1:size(AdditionalTESec,2)
        if ind_rep == 1
            for ind_sl = 1:Nslices
                n_exp = (ind_rep-1)*Nslices+ind_sl;
                AddTRTime = sum(AddTRList(1:Nslices)).*1000; % in second --> *1000
                SeqDur = SliceTime  * (Nslices-1); % nslice -1 due to saturation
                SeqDur = SeqDur + sum(EffectiveTEList(1:Nslices-1));
                CaclSliceTRtime(n_exp) = AddTRTime + SeqDur; 
            end
        else
            for ind_sl = 1:Nslices
                n_exp = (ind_rep-1)*Nslices+ind_sl;
                AddTRTime = sum(AddTRList(n_exp-(Nslices-1):n_exp)).*1000;
                SeqDur = SliceTime  * (Nslices-1); % nslice -1 due to saturation
                SeqDur = SeqDur + sum(EffectiveTEList(n_exp-(Nslices-1):n_exp-1));
                CaclSliceTRtime(n_exp) = AddTRTime + SeqDur;
            end
        end
    end
    EffecectiveTRcalc = mean(reshape(CaclSliceTRtime,Nslices,size(CaclSliceTRtime,2)/Nslices),1)./1000; % Check if constant between slices
else
    d0 = Delays(0+1);
    d6 = Delays(6+1);
    d5 = Delays(5+1);
    u15 = 0.000015;
    EffecectiveTRcalc = d0+d6+u15+AddTRList+d5+d4+d8+(p0*(1-ExcPulse1_Rpfac));
end
% figure()
% hold on
% plot(EffectiveTRList,'r')
% plot(EffecectiveTRcalc,'b')

%% Compute shapes
dt = 0.5*8e-6; % 4 useconde
Ndt_DwD4 = round(DwD4/dt);

for nShape = 1:DwNShapes
    ShapeStru.(['DWGS' num2str(nShape-1) 'R']) = ReadPV360Param(data_path_pv,['DWGS' num2str(nShape-1) 'R']);
    ShapeStru.(['DWGS' num2str(nShape-1) 'P']) = ReadPV360Param(data_path_pv,['DWGS' num2str(nShape-1) 'P']);
    ShapeStru.(['DWGS' num2str(nShape-1) 'S']) = ReadPV360Param(data_path_pv,['DWGS' num2str(nShape-1) 'S']);
    
    % Avoid error if shape is one points
    if numel(ShapeStru.(['DWGS' num2str(nShape-1) 'R'])) == 1
        Ndt = round(min(DwD4)/dt);
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
Ndt_d1 = round(d1/dt);
Ndt_d2 = round(d2/dt);
Ndt_d3 = round(d3/dt);
Ndt_d4 = round(d4/dt);
Ndt_d7 = round(d7/dt);
Ndt_d8 = round(d8/dt);
Ndt_d9 = round(d9/dt);
Ndt_d11 = round(d11/dt);
Ndt_AddTE = round(AdditionalTESec/dt);
Ndt_DwD1 = round(DwD1/dt);

Ndt_Epi20u = round(20*1e-6/dt);
Ndt_Epi10u = round(10*1e-6/dt);
Ndt_EpiD4 = round(EpiD4/dt);
Ndt_EpiD5 = round(EpiD5/dt);
Ndt_EpiD11 = round(EpiD11/dt);
% EpiEchoPosition=ReadPV360Param(data_path_pv, 'PVM_EpiEchoPosition');
% NEpiEchoCenter = round(EpiNEchoes/100*EpiEchoPosition);
DoubleSampling=ReadPV360Param(data_path_pv, 'PVM_EpiCombine') ;
NEpiEchoCenter=ReadPV360Param(data_path_pv, 'PVM_EncCentralStep1');
if strcmp(DoubleSampling,'yes')==1
   NEpiEchoCenter = (NEpiEchoCenter.*2)-1; 
end

%% Interpolate EpiShapes
list_shape = [1 3 18 14 9 11 8 2 4 13 12 7 10 5 6];
list_delays = [4 4 4 4 5 5 11 5 5 11 5 5 5 4 4];

for indshape = 1:size(list_shape,2)
    SNameTmp = ['s' num2str(list_shape(indshape))];
    DNameTmp =['Ndt_EpiD' num2str(list_delays(indshape))];
    EpiShape.(SNameTmp)= interp1(0:(eval(DNameTmp)-1)/(size(EpiShape.(SNameTmp),2)-1):(eval(DNameTmp)-1),EpiShape.(SNameTmp),0:(eval(DNameTmp)-1),'linear');%
    EpiShape.(SNameTmp)= EpiShape.(SNameTmp)'.*100;
%     EpiShape.(SNameTmp)= EpiShape.(SNameTmp)'.*0;
end

%% Check aquisition time
% CalcAcqTime = 20*1e-6 + EpiD4 +EpiD5 +EpiD11 + ...
%     (EpiD5 + EpiD11 +EpiD5 + EpiD11).*(EpiNEchoes)+....
%      EpiD5 + EpiD4; % + EpiD1
% EpiModuleTime=ReadPV360Param(data_path_pv, 'PVM_EpiModuleTime')*1e-3 ;


%% Compute main sequence part timing
Ndt_preref = Ndt_p0half + Ndt_d3 + Ndt_d1 + Ndt_d4 + Ndt_d7 +...
    Ndt_AddTE + (Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1).*DiffusionExperimentList  +...
    Ndt_d11 + Ndt_d9 + Ndt_d4 + Ndt_p1half;
Ndt_postref = Ndt_p1half + Ndt_d9 + Ndt_d4 +...
    (Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1).*DiffusionExperimentList + Ndt_AddTE +  ...
    Ndt_d11 + Ndt_d2 + ...
    Ndt_Epi20u + Ndt_Epi10u + Ndt_EpiD4 +Ndt_EpiD5 +Ndt_EpiD11 + ... % epi mod
    (Ndt_EpiD5 + Ndt_EpiD11 + Ndt_EpiD5 + Ndt_EpiD11).*((NEpiEchoCenter-1)/2); % loop
Ndt_tot = Ndt_preref + Ndt_postref;
% offset = abs(Ndt_preref - Ndt_postref)*dt*1000; %offset

%% Create sequences part for b-calculation
Seq_part(1,:) = Ndt_p0half + Ndt_d3 + Ndt_d1 + Ndt_d4 + Ndt_d7 + Ndt_AddTE;
Seq_part(2,:) = Seq_part(1,:) + (Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1).*DiffusionExperimentList;
Seq_part(3,:) = Seq_part(2,:) + Ndt_d11 + Ndt_d9 + Ndt_d4 + Ndt_p1half + Ndt_p1half + Ndt_d9 + Ndt_d4;
Seq_part(4,:) = Seq_part(3,:) + (Ndt_DwD1 + Ndt_DwD4 +Ndt_DwD1).*DiffusionExperimentList + Ndt_AddTE;
Seq_part(5,:) = Seq_part(4,:) + Ndt_d11 + Ndt_d2 + ...
    Ndt_Epi20u +  Ndt_Epi10u + Ndt_EpiD4 +Ndt_EpiD5 +Ndt_EpiD11 + ... % epi mod
    (Ndt_EpiD5 + Ndt_EpiD11 + Ndt_EpiD5 + Ndt_EpiD11).*((NEpiEchoCenter-1)/2);

xps.gwf_x = zeros(td1,max(Ndt_tot)); %Use max time step
xps.gwf_y = zeros(td1,max(Ndt_tot));
xps.gwf_z = zeros(td1,max(Ndt_tot));

%% Check echo time calculations
NdtTE = Ndt_tot*dt;
% figure()
% hold on
% plot(NdtTE,'b')
% plot(EffectiveTEList,'r')

unit2si = 1/100*PVM_GradCalConst*2*pi/gamma*1e3;
% Recompute g0
g0_exc = - g1*Ndt_d1/(Ndt_p0half+Ndt_d3/2*(d4/d3));
%% Asociate delays with gradients
G.exc_x = [0*ones(Ndt_p0half,1); 0*ones(Ndt_d3,1); 0*ones(Ndt_d1,1); 0*ones(Ndt_d4,1); 0*ones(Ndt_d7,1)]*unit2si;
G.exc_y = [0*ones(Ndt_p0half,1); 0*ones(Ndt_d3,1); 0*ones(Ndt_d1,1); 0*ones(Ndt_d4,1); 0*ones(Ndt_d7,1)]*unit2si;
G.exc_z = [g0_exc*ones(Ndt_p0half,1); g0_exc/2*(d4/d3)*ones(Ndt_d3,1); g1*ones(Ndt_d1,1); 0*ones(Ndt_d4,1); 0*ones(Ndt_d7,1)]*unit2si;

G.ref_x = [0*ones(Ndt_d11,1);0*ones(Ndt_d9,1);0*ones(Ndt_d4,1); 0*ones(Ndt_p1half,1); 0*ones(Ndt_p1half,1); 0*ones(Ndt_d9,1); 0*ones(Ndt_d4,1)]*unit2si;
G.ref_y = [0*ones(Ndt_d11,1);0*ones(Ndt_d9,1);0*ones(Ndt_d4,1); 0*ones(Ndt_p1half,1); 0*ones(Ndt_p1half,1); 0*ones(Ndt_d9,1); 0*ones(Ndt_d4,1)]*unit2si;
G.ref_z = [0*ones(Ndt_d11,1);g3*ones(Ndt_d9,1);g0*ones(Ndt_d4,1); g0*ones(Ndt_p1half,1); g0*ones(Ndt_p1half,1); g3*ones(Ndt_d9,1); g0*ones(Ndt_d4,1)]*unit2si;

G.acq_x = [0*ones(Ndt_d11,1);0*ones(Ndt_d2,1); 0*ones(Ndt_Epi20u,1);0*ones(Ndt_Epi10u,1); EpiShape.s1.*ones(Ndt_EpiD4,1);EpiShape.s9.*ones(Ndt_EpiD5,1);EpiShape.s8.*ones(Ndt_EpiD11,1)]*unit2si;
G.acq_y = [0*ones(Ndt_d11,1);0*ones(Ndt_d2,1); 0*ones(Ndt_Epi20u,1);0*ones(Ndt_Epi10u,1); (EpiShape.s3 + EpiShape.s18).*ones(Ndt_EpiD4,1);EpiShape.s11.*ones(Ndt_EpiD5,1);0*ones(Ndt_EpiD11,1)]*unit2si;
G.acq_z = [0*ones(Ndt_d11,1);0*ones(Ndt_d2,1); 0*ones(Ndt_Epi20u,1);0*ones(Ndt_Epi10u,1); -EpiShape.s14.*ones(Ndt_EpiD4,1);0*ones(Ndt_EpiD5,1);0*ones(Ndt_EpiD11,1)]*unit2si;

G.acq_x = cat(1,G.acq_x,repmat([EpiShape.s2.*ones(Ndt_EpiD5,1);EpiShape.s13.*ones(Ndt_EpiD11,1);EpiShape.s12.*ones(Ndt_EpiD5,1);EpiShape.s8.*ones(Ndt_EpiD11,1)]*unit2si,(NEpiEchoCenter-1)/2,1));
G.acq_y = cat(1,G.acq_y,repmat([EpiShape.s4.*ones(Ndt_EpiD5,1);0*ones(Ndt_EpiD11,1);EpiShape.s7.*ones(Ndt_EpiD5,1);0*ones(Ndt_EpiD11,1)]*unit2si,(NEpiEchoCenter-1)/2,1));
G.acq_z = cat(1,G.acq_z,repmat([0*ones(Ndt_EpiD5,1);0*ones(Ndt_EpiD11,1);0*ones(Ndt_EpiD5,1);0*ones(Ndt_EpiD11,1)]*unit2si,(NEpiEchoCenter-1)/2,1));

% G.acq_y = G.acq_y.*-100; % reverse acquisition

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
    
%         %% Check cumsum
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

%% invert X and Y
    x_tmp = G.x;
    G.x = G.y;
    G.y = -G.z;
    G.z = x_tmp;
%% end of invertion    
    
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

% xps.te = EffectiveTEList;
xps.te = NdtTE';

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


