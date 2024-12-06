function Spec = fPT(tauc,S,P)

omega0H = P.gammaH*P.B0;
omega0C = P.gammaC*P.B0;
omegaR = 2*pi*P.nuR;
omega1H_CP = 2*pi*P.nu1H_CP;
omega1C_rampCPvector = 2*pi*P.nu1C_rampCPvector;
Deltaomega1_CPvector = omega1H_CP-omega1C_rampCPvector;
taus = P.taus;

S2vector = S.^2;
taufvector = 1./(1./tauc + 1./P.taus);

[taufarray2,S2array2] = ndgrid(taufvector,S2vector);
%[taufarray2,Sarray2] = ndgrid(taufvector,Svector);

[taufarray3,S2array3,Deltaomega1_CParray3] = ndgrid(taufvector,S2vector,Deltaomega1_CPvector);


J0array2 = fSpecDensTwostepMAS(0,omegaR,taufarray2,taus,S2array2);
JLarmorHarray2 = fSpecDensTwostepMAS(omega0H,omegaR,taufarray2,taus,S2array2);
JLarmorCarray2 = fSpecDensTwostepMAS(omega0C,omegaR,taufarray2,taus,S2array2);
JnutH_CParray2 = fSpecDensTwostepMAS(omega1H_CP,omegaR,taufarray2,taus,S2array2);
JDeltanu1_CParray3 = fSpecDensTwostepMAS(Deltaomega1_CParray3,omegaR,taufarray3,taus,S2array3);
JDeltanu1_CParray2 = squeeze(mean(JDeltanu1_CParray3,3));

R1Carray2 = P.gammaC^2*P.bC.^2.*JLarmorCarray2;
R2Harray2 = P.gammaH^2*P.bH.^2*.5.*(J0array2 + JLarmorHarray2);
R2Carray2 = P.gammaC^2*P.bC.^2*.5.*(J0array2 + JLarmorCarray2);
R1rhoHarray2 = P.gammaH^2*P.bH.^2*.5.*(JnutH_CParray2 + JLarmorHarray2);
RCHarray2 = P.gammaC^2*P.bC.^2*.5.*JDeltanu1_CParray2;

Spec.DP = 1 - exp(-P.d1*R1Carray2);

Spec.CP = P.gammaH/P.gammaC*(exp(-P.tauCP.*R1rhoHarray2)-exp(-P.tauCP.*RCHarray2))...
    ./(1-R1rhoHarray2./RCHarray2);

Spec.INEPT = P.gammaH/P.gammaC*P.nH*sin(2*pi*P.JCH*P.tau).*sin(2*pi*P.JCH*P.tauprim)...
    .*cos(2*pi*P.JCH*P.tauprim).^(P.nH-1).*exp(-2*P.tau.*R2Harray2).*exp(-2*P.tauprim.*R2Carray2);

