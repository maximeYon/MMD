%%read time domain data and pick out relevant signal points
eval(['load ' num2str(expno(nexp)) '/NMRacqu2s'])

td = 256*ceil(NMRacqus.td/256);

tdim.z = td/2;
tdim.y = NMRacqus.l2;
tdim.x = NMRacqus.l3;
tdim.td1 = NMRacqus.l1;

tdim.tot = tdim.z*tdim.y*tdim.x*tdim.td1;

dw = 1/NMRacqus.sw_h/2;
t = 2*dw*(0:(tdim.z-1))';
t = t-t(tdim.z/2+1);

fid = fopen([num2str(expno(nexp)) '/ser'],'r','ieee-le');
Std1 = fread(fid,[2,tdim.tot],'long')';
Std1 = Std1(:,1) + i*Std1(:,2);
fclose(fid);

Std1 = reshape(Std1,tdim.z,tdim.y,tdim.x,tdim.td1);

if exist('td1start') == 0
    td1start = 1;
end

grpdlycount = (1:tdim.z)'/tdim.z - floor(tdim.z/2);
zeroshiftfun = exp(-i*((NMRacqus.grpdly-NMRacqus.de*1e-6/2/NMRacqus.dw)*2*pi*grpdlycount));
zeroshiftfun = flipdim(zeroshiftfun,1);
zeroshiftfun = fftshift(zeroshiftfun,1);
zeroshiftfun = repmat(zeroshiftfun,[1 tdim.y tdim.x tdim.td1]);

Std1 = fft(Std1,tdim.z,1);
Std1 = zeroshiftfun.*Std1;
Std1 = ifft(Std1,tdim.z,1);

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

g.z = NMRacqus.cnst13*Gmax/100;
k.z = gamma*g.z*t/2/pi;

g.x = linspace(-1,1,tdim.x+1)'*Gmax*NMRacqus.cnst21/100; g.x = g.x(1:tdim.x);
g.y = linspace(-1,1,tdim.y+1)'*Gmax*NMRacqus.cnst22/100; g.y = g.y(1:tdim.y);

t_phasenc = NMRacqus.d33 + NMRacqus.d32;
k.x = gamma*g.x'*t_phasenc/2/pi;
k.y = gamma*g.y'*t_phasenc/2/pi;

%%smoothing
[karray.z,karray.y,karray.x,dummy] = ndgrid(k.z,k.y,k.x,1:tdim.td1);
lbfun = exp(-lb.^2*pi.^2*(karray.z.^2 + karray.y.^2 + karray.x.^2));

Std1 = Std1 - mean(mean(mean(mean(Std1([1:round(.1*tdim.z) round(.9*tdim.z):tdim.z],:,:,:)))));
Std1 = lbfun.*Std1;


if exist('six') == 0
    nudim.x = NMRacqus.l3;
else
    nudim.x = six;
end
if exist('siy') == 0
    nudim.y = NMRacqus.l2;
else
    nudim.y = siy;
end
if exist('si') == 0
    nudim.z = NMRacqus.td/2;
else
    nudim.z = si;
end

if nudim.z < tdim.z
    index = .5*(tdim.z-nudim.z) + (1:nudim.z);
    Std1 = Std1(index,:,:,:);
    karray.z = karray.z(index,:,:,:);
    karray.y = karray.y(index,:,:,:);
    karray.x = karray.x(index,:,:,:);
end


%%Fourier transform
Stemp = Std1;
if nudim.z > tdim.z
    Nzeros = nudim.z - tdim.z;
    Stemp2 = zeros(nudim.z,tdim.y,tdim.x,tdim.td1);
    Stemp2(Nzeros/2+(1:tdim.z),:,:,:) = Std1;
    Stemp = Stemp2;
end
Stemp = fftshift(fft(ifftshift(Stemp,1),[],1),1);
if nudim.y > tdim.y
    Nzeros = nudim.y - tdim.y;
    Stemp2 = zeros(nudim.z,nudim.y,tdim.x,tdim.td1);
    Stemp2(:,Nzeros/2+(1:tdim.y),:,:) = Stemp;
    Stemp = Stemp2;
end
Stemp = fftshift(fft(ifftshift(Stemp,2),[],2),2);
if nudim.x > tdim.x
    Nzeros = nudim.x - tdim.x;
    Stemp2 = zeros(nudim.z,nudim.y,nudim.x,tdim.td1);
    Stemp2(:,:,Nzeros/2+(1:tdim.x),:) = Stemp;
    Stemp = Stemp2;
end
Itd1 = fftshift(fft(ifftshift(Stemp,3),[],3),3);
clear Stemp Stemp2

Imax = abs(Itd1);
Imax(nudim.z/2+(-2:4),nudim.y/2+(0:2),nudim.x/2+(0:2),:) = 0;
Imax = max(reshape(Imax,numel(Imax),1));

r.x = .5/(k.x(2)-k.x(1))*linspace(-1,1,nudim.x+1)'; r.x = r.x(1:nudim.x);
r.y = .5/(k.y(2)-k.y(1))*linspace(-1,1,nudim.y+1)'; r.y = r.y(1:nudim.y);
r.z = .5/(k.z(2)-k.z(1))*linspace(-1,1,nudim.z+1)'; r.z = r.z(1:nudim.z);

resolution.x = abs(r.x(2)-r.x(1));
resolution.y = abs(r.y(2)-r.y(1));
resolution.z = abs(r.z(2)-r.z(1));

FOV.x = max(r.x) - min(r.x);
FOV.y = max(r.y) - min(r.y);
FOV.z = max(r.z) - min(r.z);
