clear all
%clc

B1ref1H = 187/2;             % kHz
plwref1H = 60;                     % Watt
B1ref13C = 235/2;           % kHz
plref13C = 200;                   % Watt

spinrate = 12e3;

B1hp1H = 100;   % kHz
B1dec = 80;                   % kHz
B1hp13C = B1hp1H;   % kHz
B1cp13C = 75;          % kHz

p90hp = 1000/B1hp1H/4                    % microsecs
p180dec = 1000/B1dec/2; %us
pcpd2 = p180dec
%pl180dec = 7.35;                 % microsec
%B1dec = 1/(2*pl180dec*0.001)   % kHz

plw1 = plref13C*(B1hp13C/B1ref13C)^2   % Watt
plw21 = plref13C*(B1cp13C/B1ref13C)^2  % Watt

plw2 = plwref1H*(B1hp1H/B1ref1H)^2                % Watt
plw12 = plwref1H*(B1dec/B1ref1H)^2                % Watt


return
B1hp13C = 1/(4*pl90hp*0.001);   % kHz
B1cp13C = 0.9*B1hp13C;          % kHz
B1spinlock13C = 30;   % kHz

%pl180rec = 1/(spinrate*9*2)*10^6 % microsecs
acqtime = (pl180dec/20.4)^2*0.4
% B1cp = 63.5;

% B1rec = 45;
B1rec = 1/(2*pl180rec*0.001);   % kHz

PL2dB = plrefdB + 20*log10(B1ref1H/B1hp1H);    % dB
PL12dB = plrefdB + 20*log10(B1ref1H/B1dec);    % dB
PL22dB = PL2dB;
PL32dB = plrefdB + 20*log10(B1ref1H/B1rec);    % dB
SPW0 = PLW2
PLW32 = plref*(B1rec/B1ref1H)^2                % Watt


% B1hp13C = 80;
B1hp13C = 1/(4*pl90hp*0.001);   % kHz
B1cp13C = 0.9*B1hp13C;          % kHz
B1spinlock13C = 30;   % kHz
plref2 = 150;                   % Watt
plref2dB = -21.76;                   % dB

PL1dB = plref2dB + 20*log10(B1ref13C/B1hp13C);   % dB
PL21dB = plref2dB + 20*log10(B1ref13C/B1cp13C);  % dB
PLW1 = plref2*(B1hp13C/B1ref13C)^2   % Watt
PLW21 = plref2*(B1cp13C/B1ref13C)^2  % Watt
PLW13 = plref2*(B1spinlock13C/B1ref13C)^2  % Watt
