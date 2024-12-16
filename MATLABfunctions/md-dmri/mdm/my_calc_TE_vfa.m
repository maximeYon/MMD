function [TEeq] = my_calc_TE_vfa(etl,prepEsp,esp,flipangle)
%% Calculation of TE for massively multidimensionnal diffusion
% Maxime Yon 2022

% etl=40; % Number of echoes (acquired + discarded)
% prepEsp = 0.01; % diffusion spin echo duration (s)
% esp = 0.002; %echo spacing (s)
T2 = 0.030; % s
T1 = 0.900; % s
% flipangle = [2.20591678769962,1.63580219786576,1.40997184351612,1.05685133748620,0.791962980880645,0.698131700797732,0.698131700797732,0.698131700797732,0.698131700797732,0.698131700797732,0.698131700797732,0.708215928402292,0.737755614511753,0.766875934666186,0.799673534137448,0.833915850374733,0.876421968099477,0.926497418699185,0.990138667524560,1.06624673139625,1.16234544765635,1.25844416391645,1.35454288017654,1.45064159643664,1.54674031269673,1.64283902895683,1.73893774521693,1.83503646147702,1.93113517773712,2.02723389399721,2.12333261025731,2.21943132651741,2.31553004277750,2.41162875903760,2.50772747529770,2.60382619155779,2.69992490781789,2.79602362407798,2.89212234033808,2.98822105659818];
% % in radian

%% Intitialization
P = zeros(3,2*etl+2);		% Allocate all known states, 2 per echo.
P(3,1)=1;			% Initial condition/equilibrium.
Pcoh = zeros(3,2*etl+2);		% Allocate all known states, 2 per echo.
Pcoh(3,1)=1;			% Initial condition/equilibrium.
s = zeros(1,etl+1);		% Allocate signal vector to store.
fcoh = zeros(1,etl+1);		% Allocate signal vector to store.
echostart = 2; % with spin echo preparation

%% -- Preparation module for diffusion echo
P = epg_rf(P,pi/2,pi/2);	% Do 90 tip.
P = epg_grelax(P,T1,T2,prepEsp/2,1,0,1,0);
P = epg_rf(P,pi,pi);   % -- Refoc. RF
P = epg_grelax(P,T1,T2,prepEsp/2,1,0,1,0);
s(1) = abs(P(1,1));  	% Signal is F0 state.

Pcoh = epg_rf(Pcoh,pi/2,pi/2);	% Do 90 tip.
Pcoh = epg_grelax(Pcoh,Inf,Inf,prepEsp/2,1,0,1,0);
Pcoh = epg_rf(Pcoh,pi,pi);   % -- Refoc. RF
Pcoh = epg_grelax(Pcoh,Inf,Inf,prepEsp/2,1,0,1,0);

%% VFA echo train
% Pcoh = P;
faind = 1;
for ech=echostart:etl+1
    P = epg_grelax(P,T1,T2,esp/2,1,0,1,0);   % -- Left crusher
    Pcoh = epg_grelax(Pcoh,Inf,Inf,esp/2,1,0,1,1);   % -- Left crusher
    P = epg_rf(P,abs(flipangle(faind)),pi);   % -- Refoc. RF
    Pcoh = epg_rf(Pcoh,abs(flipangle(faind)),pi);   % -- Refoc. RF
    P = epg_grelax(P,T1,T2,esp/2,1,0,1,0);   % -- Right crusher
    Pcoh = epg_grelax(Pcoh,Inf,Inf,esp/2,1,0,1,1);   % -- Right crusher
    
    s(ech) = abs(P(1,1));  	% Signal is F0 state.
    fcoh(ech) = abs(Pcoh(1,1));  	% Signal is F0 state.
    faind = faind + 1;
end
% figure(1)
% plot(-T2*log(s./fcoh))
% drawnow

TEeq = -T2*log(s./fcoh);
TEeq = TEeq(2:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    Function needeed                                    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FpFmZ,RR] = epg_rf(FpFmZ,alpha,phi)
%
%	Propagate EPG states through an RF rotation of 
%	alpha, with phase phi (both radians).
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%

% -- From Weigel at al, JMR 205(2010)276-285, Eq. 8.

RR = [(cos(alpha/2))^2 exp(2*1i*phi)*(sin(alpha/2))^2 -1i*exp(1i*phi)*sin(alpha);
      exp(-2*1i*phi)*(sin(alpha/2))^2 (cos(alpha/2))^2 1i*exp(-1i*phi)*sin(alpha);
      -1i/2*exp(-1i*phi)*sin(alpha) 1i/2*exp(1i*phi)*sin(alpha)      cos(alpha)];


FpFmZ = RR * FpFmZ;
end

%function [FpFmZ,EE,BV] = epg_grelax(FpFmZ,T1,T2,T,kg,D,Gon,noadd)
%
%	Propagate EPG states through a period of relaxation, and
%	diffusion over an interval T, with or without a gradient.
%	Leave last 3 blank to exclude diffusion effects.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		T1,T2 = Relaxation times (s)
%		T = Time interval (s)
%		kg = k-space traversal due to gradient (rad/m)
%		D = Diffusion coefficient (m^2/s)
%		Gon = 0 if no gradient on, 1 if gradient on
%			(gradient will advance states at the end.)
%		noadd=1 to not add higher-order states - see epg_grad.m
%
%	OUTPUT:
%		FpFmZ = updated F+, F- and Z states.
%		EE = decay matrix, 3x3 = diag([E2 E2 E1]);
%		BV = b-value matrix, 3xN (see FpFmZ) of attenuations.

function [FpFmZ,EE,BV] = epg_grelax(FpFmZ,T1,T2,T,kg,D,Gon,noadd)

if (nargin < 8); noadd=0; end	% Default is to add states.

E2 = exp(-T/T2);
E1 = exp(-T/T1);

EE = diag([E2 E2 E1]);		% Decay of states due to relaxation alone.
RR = 1-E1;			% Mz Recovery, affects only Z0 state, as 
				% recovered magnetization is not dephased.


FpFmZ = EE * FpFmZ;		% Apply Relaxation
FpFmZ(3,1) = FpFmZ(3,1)+RR;	% Recovery  ( here applied before diffusion,
				% but could be after or split.)
%E1
%E2
%EE
%RR

%disp('(In epg_grelax): Before diffusion:');
%FpFmZ

if (nargin > 4)			% Model Diffusion Effects

  Findex = 0:length(FpFmZ(1,:))-1;	% index of states, 0...N-1
  bvalZ = ((Findex)*kg).^2*T;		% diffusion  for Z states, assumes that
					% the Z-state has to be refocused, so
					% this models "time between gradients"
	
	% For F states, the following models the additional diffusion time
	% (Findex) and the fact that the state will change if the gradient is
	% on (0.5*Gon), then the additional diffusion *during* the gradient, 
	% ... Gon*kg^2/12 term.

  bvalp = ((( Findex+.5*Gon)*kg).^2+Gon*kg^2/12)*T;	% for F+ states
  bvalm = (((-Findex+.5*Gon)*kg).^2+Gon*kg^2/12)*T;	% for F- states

		

  FpFmZ(1,:) = FpFmZ(1,:) .* exp(-bvalp*D);	% diffusion on F+ states
  FpFmZ(2,:) = FpFmZ(2,:) .* exp(-bvalm*D);	% diffusion on F- states
  FpFmZ(3,:) = FpFmZ(3,:) .* exp(-bvalZ*D);	% diffusion of Z states.
 
  BV = [bvalp; bvalm; bvalZ];	% For output. 
end

%disp('(In epg_grelax): After diffusion:');
%FpFmZ

if (Gon==1)
  FpFmZ = epg_grad(FpFmZ,noadd);	% Advance states.
end

%disp('(In epg_grelax): After advancement of states:');
%FpFmZ
end

function [FpFmZ] = epg_grad(FpFmZ,noadd)
%
%	Propagate EPG states through a "unit" gradient.
%	
%	INPUT:
%		FpFmZ = 3xN vector of F+, F- and Z states.
%		noadd = 1 to NOT add any higher-order states - assume
%			that they just go to zero.  Be careful - this
%			speeds up simulations, but may compromise accuracy!
%

if (nargin < 2); noadd=0; end	% Add by default.  

% Gradient does not affect the Z states.

if (noadd==0)
  FpFmZ = [FpFmZ [0;0;0]];	% Add higher dephased state.
end

FpFmZ(1,:) = circshift(FpFmZ(1,:),[0 1]);	% Shift Fp states.
FpFmZ(2,:) = circshift(FpFmZ(2,:),[0 -1]);	% Shift Fm states.
FpFmZ(2,end)=0;					% Zero highest Fm state.
FpFmZ(1,1) = conj(FpFmZ(2,1));			% Fill in lowest Fp state.
end
