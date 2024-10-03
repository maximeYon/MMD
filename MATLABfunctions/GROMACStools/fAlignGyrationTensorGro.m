function gro_out = fAlignGyrationTensorGro(gro_in)

gro = gro_in;

com.x = mean(gro.x);
com.y = mean(gro.y);
com.z = mean(gro.z);

com_neg.x = -com.x;
com_neg.y = -com.y;
com_neg.z = -com.z;

GyrationTensor = fGyrationTensorGro(gro);

euler.alpha = 0;
euler.beta = -acos(GyrationTensor.lambda33vec(3));
euler.gamma = -atan2(GyrationTensor.lambda33vec(2),GyrationTensor.lambda33vec(1));
gro = fTranslateGro(fRotateGro(fTranslateGro(gro,com_neg),euler),com);

euler.alpha = 0;
euler.beta = 0;
euler.gamma = -atan2(GyrationTensor.lambda22vec(2),GyrationTensor.lambda22vec(1));
%euler.gamma = -pi-atan2(GyrationTensor.lambda11vec(2),GyrationTensor.lambda11vec(1));
gro = fTranslateGro(fRotateGro(fTranslateGro(gro,com_neg),euler),com);

gro_out = gro;