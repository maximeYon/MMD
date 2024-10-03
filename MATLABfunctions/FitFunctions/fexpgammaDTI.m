function Y = fexpDTI(Pin,Xin,Pnorm,Xnorm,Ynorm,Gnorm_LF);Pin = Pin.*Pnorm;Xin = Xin*Xnorm;Y0 = Pin(1);lambda.x = Pin(2);lambda.y = Pin(3);lambda.z = Pin(4);A.alpha = Pin(5);A.beta = Pin(6);A.gamma = Pin(7);sigma.x = Pin(8);sigma.y = Pin(9);sigma.z = Pin(10);R.gamma = [    cos(A.gamma) sin(A.gamma) 0    -sin(A.gamma) cos(A.gamma) 0    0 0 1];R.beta = [    cos(A.beta) 0 -sin(A.beta)    0 1 0    sin(A.beta) 0 cos(A.beta)];R.alpha = [    cos(A.alpha) sin(A.alpha) 0    -sin(A.alpha) cos(A.alpha) 0    0 0 1];R.mat = R.gamma*R.beta*R.alpha;R.inv = inv(R.mat);Gnorm_PAS.x = R.inv(1,1)*Gnorm_LF.x + R.inv(1,2)*Gnorm_LF.y + R.inv(1,3)*Gnorm_LF.z;Gnorm_PAS.y = R.inv(2,1)*Gnorm_LF.x + R.inv(2,2)*Gnorm_LF.y + R.inv(2,3)*Gnorm_LF.z;Gnorm_PAS.z = R.inv(3,1)*Gnorm_LF.x + R.inv(3,2)*Gnorm_LF.y + R.inv(3,3)*Gnorm_LF.z;b_PAS.x = Gnorm_PAS.x.^2.*Xin;b_PAS.y = Gnorm_PAS.y.^2.*Xin;b_PAS.z = Gnorm_PAS.z.^2.*Xin;Y = Y0.*((1 + b_PAS.x.*(sigma.x.^2)./lambda.x).^(-1*(lambda.x./sigma.x).^2))...    .*((1 + b_PAS.y.*(sigma.y.^2)./lambda.y).^(-1*(lambda.y./sigma.y).^2))...    .*((1 + b_PAS.z.*(sigma.z.^2)./lambda.z).^(-1*(lambda.z./sigma.z).^2));Y = Y/Ynorm;