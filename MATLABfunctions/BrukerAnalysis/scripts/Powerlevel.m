clear all

plref = 0.79;
B1ref = 80e3;

pl = 2.268;
B1 = B1ref/10^((pl - plref)/20)

pl = plref + 20*log10(B1ref/B1)
