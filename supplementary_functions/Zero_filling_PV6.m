function [kFieldZF] = Zero_filling_PV6(kField,ZF)

%% zero filling
ZF = max([ZF size(kField,1)]);
kFieldZF = zeros(ZF, size(kField,2));
trunc = ZF-size(kField,1);
if mod(trunc,2)==1
    kFieldZF(round(trunc/2):end-round(trunc/2),:)=kField;
else
    kFieldZF(round(trunc/2)+1:end-round(trunc/2),:)=kField;
end

















