eqindx = strfind(eq.x,'.x');

eq.y = eq.x;
eq.z = eq.x;

for neq = 1:length(eqindx)
    eq.y(eqindx(neq)+[0 1]) = '.y';
    eq.z(eqindx(neq)+[0 1]) = '.z';
end

eval(eq.x)
eval(eq.y)
eval(eq.z)
