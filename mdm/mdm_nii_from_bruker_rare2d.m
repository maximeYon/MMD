function res = mdm_nii_from_bruker_rare2d(proc,path)
%%read time domain data and pick out relevant signal points

eval(['load ' path '/NMRacqu2s' ])
eval(['load ' path '/NMRacqus' ])

td = 256*ceil(NMRacqus.td/256);
td1 = NMRacqu2s.td;
dw = 1/NMRacqus.sw_h/2;
t = 2*dw*(0:(td/2-1))';

tdim.tot = td/2*td1;

fid = fopen([path '/ser'],'r','ieee-le');
Std1 = fread(fid,[2,tdim.tot],'long')';
Std1 = Std1(:,1) + i*Std1(:,2);
fclose(fid);

Std1 = reshape(Std1,td/2,td1);
S = Std1(:,1);

grpdlycount = (1:(td/2))'/(td/2) - floor(td/2/2);
zeroshiftfun = exp(-i*(NMRacqus.grpdly*2*pi*grpdlycount));
zeroshiftfun = flipdim(zeroshiftfun,1);
zeroshiftfun = fftshift(zeroshiftfun,1);
S = fft(S,td/2,1);
S = zeroshiftfun.*S;
S = ifft(S,td/2,1);
%figure(1), clf, plot(1:td/2,real(S),1:td/2,imag(S)), return

tdim.i = NMRacqus.d40/NMRacqus.dw/2;
tdim.j = NMRacqus.l1;
tdim.ramp = NMRacqus.d32/dw/2;
tdim.fid = 2+(NMRacqus.d32+NMRacqus.d41)/dw/2;
tdim.phase = 2*NMRacqus.d32/dw/2;
tdim.slice = 4+NMRacqus.p12/1e6/dw/2;
tdim.read = NMRacqus.d40/dw/2;
tdim.init = tdim.fid + 2*tdim.phase + tdim.slice;
tdim.repeat = 2*tdim.phase + tdim.slice + tdim.read;

tdim.count.i = 1:tdim.i;
tdim.count.j = 0:(tdim.j-1);

[tdim.count.array.i,tdim.count.array.j] = ndgrid(tdim.count.i,tdim.count.j);
index.S = tdim.init + tdim.count.array.j*tdim.repeat + tdim.count.array.i;

S = reshape(S(index.S),tdim.i,tdim.j);
%figure(1), clf, plot(1:tdim.i,real(S(:,(tdim.j/2+1))),1:tdim.i,imag(S(:,(tdim.j/2+1)))), return


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

%%read i
g.i = sqrt(NMRacqus.cnst11^2+NMRacqus.cnst12^2+NMRacqus.cnst13^2)*Gmax/100;
t_read = 2*dw*(0:(tdim.i-1))'; t_read = t_read-t_read(tdim.i/2+1);
k.i = gamma*g.i*t_read/2/pi;

%%phase j
g.j = linspace(-1,1,tdim.j+1)'*sqrt(NMRacqus.cnst21^2+NMRacqus.cnst22^2+NMRacqus.cnst23^2)*Gmax/100;
g.j = g.j(1:tdim.j);
t_phasenc = 0*NMRacqus.d33 + NMRacqus.d32;
k.j = gamma*g.j*t_phasenc/2/pi;

%%smoothing
lbfun.i = exp(-(proc.lb*pi*k.i).^2);
lbfun.j = exp(-(proc.lb*pi*k.j).^2);
%figure(1), clf, plot(k.i,lbfun.i,'o',k.j,lbfun.j,'s'), return

S = S - mean(mean(S([1:round(.1*tdim.i) round(.9*tdim.i):tdim.i],:)));

Stemp = S(:,tdim.j/2+1);
%figure(1), clf, plot(1:tdim.i,real(Stemp),1:tdim.i,imag(Stemp)), return
Stemp = lbfun.i.*Stemp;
%figure(1), clf, plot(k.i,real(Stemp),k.i,lbfun.i*max(real(Stemp))), return


if isfield(proc,'si') == 0
    nudim.i = tdim.i;
else
    nudim.i = proc.si;
end
if isfield(proc,'si1') == 0
    nudim.j = tdim.j;
else
    nudim.j = proc.si1;
end

index.tdim.i = 1:tdim.i;
if nudim.i<tdim.i
    index.tdim.i = .5*(tdim.i-nudim.i) + (1:nudim.i);
    Stemp = Stemp(index.tdim.i);
    k.i = k.i(index.tdim.i);
    lbfun.i = lbfun.i(index.tdim.i);
    %figure(1), clf,plot(k.i,real(Stemp),k.i,lbfun.i*max(real(Stemp))), return
end
index.tdim.j = 1:tdim.j;
if nudim.j<tdim.j
    index.tdim.j = .5*(tdim.j-nudim.j) + (1:nudim.j);
    k.j = k.j(index.tdim.j);
    lbfun.j = lbfun.j(index.tdim.j);
end

Itemp = fftshift(fft(Stemp,nudim.i));
%figure(1), clf, plot(abs(Itemp)), return

[k.array.i,k.array.j] = ndgrid(k.i,k.j);
lbfun.array = exp(-proc.lb.^2*pi.^2*(k.array.i.^2 + k.array.j.^2));

S = S(index.tdim.i,index.tdim.j);

S = S.*lbfun.array;
%figure(1), clf, plot(real(squeeze(S(:,17)))), return

%%Fourier transform
if nudim.i > tdim.i
    Nzeros.i = nudim.i - tdim.i;
    S = [zeros(Nzeros.i/2,tdim.j); S; zeros(Nzeros.i/2,tdim.j)];
end
I = fftshift(fft(ifftshift(S,1),[],1),1);

if isfield(NMRacqus,'fq1')
    %fqcycle = 3480;
    %fqcycle = 1e10;
    if NMRacqus.dw == 5e-6, fqcycle = -2425;
    elseif NMRacqus.dw == 10e-6, fqcycle = -2340;
    end
    Nrep = td1/numel(NMRacqus.fq1);
    if Nrep>1, fq1list = repmat(NMRacqus.fq1,[Nrep 1]);
    else fq1list = NMRacqus.fq1(1:td1);
    end
    phcorrfun = repmat(exp(i*2*pi*fq1list(1)/fqcycle*((1:tdim.j)-tdim.j/2)),[nudim.i 1]);
    I = phcorrfun.*I;
end

if nudim.j > tdim.j
    Nzeros.j = nudim.j - tdim.j;
    I = [zeros(nudim.i,Nzeros.j/2) I zeros(nudim.i,Nzeros.j/2)];
end
I = fftshift(fft(ifftshift(I,2),[],2),2);
%figure(1), clf, imagesc(abs(I)'), set(gca,'YDir','normal'), return
r.i = .5/(k.i(2)-k.i(1))*linspace(-1,1,nudim.i+1)'; r.i = r.i(1:nudim.i);
r.j = .5/(k.j(2)-k.j(1))*linspace(-1,1,nudim.j+1)'; r.j = r.j(1:nudim.j);

resolution.i = r.i(2) - r.i(1);
resolution.j = r.j(2) - r.j(1);

Itd1 = zeros(nudim.i,nudim.j,1,td1);
for ntd1 = 1:td1
    S = Std1(:,ntd1);
    S = fft(S,td/2,1);
    S = zeroshiftfun.*S;
    S = ifft(S,td/2,1);
    S = reshape(S(index.S),tdim.i,tdim.j);
    S = S(index.tdim.i,index.tdim.j);
    S = S.*lbfun.array;
    if nudim.i > tdim.i
        S = [zeros(Nzeros.i/2,tdim.j); S; zeros(Nzeros.i/2,tdim.j)];
    end
    I = fftshift(fft(ifftshift(S,1),[],1),1);
    if isfield(NMRacqus,'fq1')
        phcorrfun = repmat(exp(i*2*pi*fq1list(ntd1)/fqcycle*((1:tdim.j)-tdim.j/2)),[nudim.i 1]);
        I = phcorrfun.*I;
    end
    if nudim.j > tdim.j
        I = [zeros(nudim.i,Nzeros.j/2) I zeros(nudim.i,Nzeros.j/2)];
    end
    I = fftshift(fft(ifftshift(I,2),[],2),2);
    Itd1(:,:,1,ntd1) = I;
end

%figure(1), clf, imagesc(abs(squeeze(Itd1(:,:,1,1)))), return
mdm_mat2nii(Itd1,[path '/in.nii'],[resolution.i resolution.j 1 1])

res = 1;


