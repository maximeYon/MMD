function [r,g,b] = fphase2rgb(phase)

r = zeros(size(phase));
g = r;
b = r;

phase = phase/2/pi;
pos = find(phase<0);
phase(pos) = phase(pos)+1;

pos = find(phase>=0 & phase<1/6);
r(pos) = 1;
g(pos) = phase(pos)*6;

pos = find(phase>=1/6 & phase<1/3);
g(pos) = 1;
r(pos) = 1-(phase(pos)-1/6)*6;

pos = find(phase>=1/3 & phase<1/2);
g(pos) = 1;
b(pos) = (phase(pos)-1/3)*6;

pos = find(phase>=1/2 & phase<2/3);
b(pos) = 1;
g(pos) = 1-(phase(pos)-1/2)*6;

pos = find(phase>=2/3 & phase<5/6);
b(pos) = 1;
r(pos) = (phase(pos)-2/3)*6;

pos = find(phase>=5/6 & phase<1);
r(pos) = 1;
b(pos) = 1-(phase(pos)-5/6)*6;
