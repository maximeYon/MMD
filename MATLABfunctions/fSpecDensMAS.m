function J = fSpecDensMAS(omega,omegaR,tauc)

J = 1/3*(fLorentzSpecDens(omega-omegaR,tauc)+fLorentzSpecDens(omega+omegaR,tauc)) + ...
    1/6*(fLorentzSpecDens(omega-2*omegaR,tauc)+fLorentzSpecDens(omega+2*omegaR,tauc));
