function dataOut=fineshift(dataIn,step)

datasize=size(dataIn);
if length(step)~=length(datasize);
    error('shift step need the same size of the input data');
end

FFTdataIn=fftshift(fft2(dataIn));


% aa = ones(datasize(1),1)*exp(1i*2*pi*step(2)*(0:datasize(2)-1)/datasize(2));
%   bb = aa.*exp(1i*2*pi*step(1)*(0:datasize(1)-1).'/datasize(1)*ones(1,datasize(2)));

dataOut=ifft2(ifftshift(FFTdataIn.*(ones(datasize(1),1)*exp(1i*2*pi*step(2)*(0:datasize(2)-1)/datasize(2)))...
    .*exp(1i*2*pi*step(1)*(0:datasize(1)-1).'/datasize(1)*ones(1,datasize(2)))));