clear all

B1ref1H = 96.9/2;             % kHz
% B1ref1H = 1/(4*10.81e-6*1000)
B1ref13C = 144.3/2;           % kHz
pl90hp = 3.1;                    % microsecs
spinrate = 5e3;

pl180rec = 1/(spinrate*9*2)*10^6 % microsecs
pl180dec = 7.35;                 % microsec
acqtime = (pl180dec/20.4)^2*0.4
% B1hp1H = 80;
B1hp1H = 1/(4*pl90hp*0.001);    % kHz
% B1cp = 63.5;
% B1dec = 25;                   % 25 kHz
B1dec = 1/(2*pl180dec*0.001);   % kHz
% B1rec = 45;
B1rec = 1/(2*pl180rec*0.001);   % kHz
plref = 33;                     % Watt
plrefdB = -15.19;                    % dB

PL2dB = plrefdB + 20*log10(B1ref1H/B1hp1H);    % dB
PL12dB = plrefdB + 20*log10(B1ref1H/B1dec);    % dB
PL22dB = PL2dB;
PL32dB = plrefdB + 20*log10(B1ref1H/B1rec);    % dB
PLW2 = plref*(B1hp1H/B1ref1H)^2                % Watt
PLW12 = plref*(B1dec/B1ref1H)^2                % Watt
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
