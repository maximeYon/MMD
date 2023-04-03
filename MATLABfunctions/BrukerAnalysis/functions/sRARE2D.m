%%read time domain data and pick out relevant signal points
td = 256*ceil(NMRacqus.td/256);
dw = 1/NMRacqus.sw_h/2;
t = 2*dw*(0:(td/2-1))';

fid = fopen([num2str(expno(nexp)) '/fid'],'r','ieee-le');
%S = fread(fid,[2,td/2],'long')';
S = fread(fid,[2,td/2],'float64')';
S = S(:,1) + 1i*S(:,2);
fclose(fid);
%figure(1), clf, plot(1:td/2,real(S),1:td/2,imag(S)), return

grpdlycount = (1:(td/2))'/(td/2) - floor(td/2/2);
zeroshiftfun = exp(-i*((NMRacqus.grpdly-.1)*2*pi*grpdlycount))*exp(-i*pi/2);
zeroshiftfun = flipdim(zeroshiftfun,1);
zeroshiftfun = fftshift(zeroshiftfun,1);
S = fft(S,td/2,1);
S = zeroshiftfun.*S;
S = ifft(S,td/2,1);
% figure(1), clf, plot(1:td/2,real(S),1:td/2,imag(S)), return

tdim.i = NMRacqus.d40/NMRacqus.dw/2;
tdim.j = NMRacqus.l1;
tdim.ramp = NMRacqus.d32/dw/2;
tdim.fid = 2+(NMRacqus.d32+NMRacqus.d41)/dw/2;
tdim.phase = 2*NMRacqus.d32/dw/2;
tdim.slice = 4+NMRacqus.p12/1e6/dw/2;
tdim.read = NMRacqus.d40/dw/2;
tdim.init = tdim.fid + 2*tdim.phase + tdim.slice + 2;
tdim.repeat = 2*tdim.phase + tdim.slice + tdim.read;

tdim.count.i = 1:tdim.i;
tdim.count.j = 0:(tdim.j-1);

[tdim.count.array.i,tdim.count.array.j] = ndgrid(tdim.count.i,tdim.count.j);
index.S = tdim.init + tdim.count.array.j*tdim.repeat + tdim.count.array.i;

S = reshape(S(index.S),tdim.i,tdim.j);
% figure(1), clf, plot(1:tdim.i,real(S(:,(tdim.j/2+1))),1:tdim.i,imag(S(:,(tdim.j/2+1)))), return


%%calculate k-space
gamma = 26.75e7;
if strcmp(NMRacqus.nuc1,'2H') == 1
    gamma = 4.1065e7;
elseif strcmp(NMRacqus.nuc1,'23Na') == 1
    gamma = 7.0761e7;
end

Gmax = 3;
if any(strcmp(NMRacqus.probhd,{'5 mm BBO BB-1H/D XYZ-GRD Z107255/0001',...
        '5 mm TXI 1H/D-13C/15N XYZ-GRD Z8588/0006'})) == 1
    Gmax = 0.5;
end

%%read
g.i = sqrt(NMRacqus.cnst11^2+NMRacqus.cnst12^2+NMRacqus.cnst13^2)*Gmax/100;
t_read = 2*dw*(0:(tdim.i-1))'; t_read = t_read-t_read(tdim.i/2+1);
k.i = gamma*g.i*t_read/2/pi;

%%phase
g.j = linspace(-1,1,tdim.j+1)'*sqrt(NMRacqus.cnst21^2+NMRacqus.cnst22^2+NMRacqus.cnst23^2)*Gmax/100;
g.j = g.j(1:tdim.j);
t_phasenc = 0*NMRacqus.d33 + NMRacqus.d32;
k.j = gamma*g.j*t_phasenc/2/pi;


%%smoothing
lbfun.i = exp(-(lb*pi*k.i).^2);
lbfun.j = exp(-(lb*pi*k.j).^2);
%figure(1), clf, plot(k.i,lbfun.i,'o',k.j,lbfun.j,'s'), return

S = S - mean(mean(S([1:round(.1*tdim.i) round(.9*tdim.i):tdim.i],:)));

Stemp = S(:,tdim.j/2+1);
%figure(1), clf, plot(1:tdim.i,real(Stemp),1:tdim.i,imag(Stemp)), return
Stemp = lbfun.i.*Stemp;
%figure(1), clf, plot(k.i,real(Stemp),k.i,lbfun.i*max(real(Stemp))), return

if exist('si') == 0
    nudim.i = tdim.i;
else
    nudim.i = si;
end
if exist('si1') == 0
    nudim.j = tdim.j;
else
    nudim.j = si1;
end

index.tdim.i = 1:tdim.i;
if nudim.i<tdim.i
    index.tdim.i = .5*(tdim.i-nudim.i) + (1:nudim.i);
    Stemp = Stemp(index.tdim.i);
    k.i = k.i(index.tdim.i);
    lbfun.i = lbfun.i(index.tdim.i);
    %figure(1), clf,plot(k.i,real(Stemp),k.i,lbfun.i*max(real(Stemp))), return
end

Itemp = fftshift(fft(Stemp,nudim.i));
%figure(1), clf, plot(abs(Itemp)), return

[k.array.i,k.array.j] = ndgrid(k.i,k.j);
lbfun.array = exp(-lb.^2*pi.^2*(k.array.i.^2 + k.array.j.^2));

S = S(index.tdim.i,:);

S = S.*lbfun.array;
%figure(1), clf, plot(real(squeeze(S(:,17)))), return

%%Fourier transform
if nudim.i > tdim.i
    Nzeros.i = nudim.i - tdim.i;
    S = [zeros(Nzeros.i/2,tdim.j); S; zeros(Nzeros.i/2,tdim.j)];
end
I = fftshift(fft(ifftshift(S,1),[],1),1);
if nudim.j > tdim.j
    Nzeros.j = nudim.j - tdim.j;
    I = [zeros(nudim.i,Nzeros.j/2) I zeros(nudim.i,Nzeros.j/2)];
end
I = fftshift(fft(ifftshift(I,2),[],2),2);

r.i = .5/(k.i(2)-k.i(1))*linspace(-1,1,nudim.i+1)'; r.i = r.i(1:nudim.i);
r.j = .5/(k.j(2)-k.j(1))*linspace(-1,1,nudim.j+1)'; r.j = r.j(1:nudim.j);

resolution.i = r.i(2) - r.i(1);
resolution.j = r.j(2) - r.j(1);
