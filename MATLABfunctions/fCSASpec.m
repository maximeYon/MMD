function Spec = fCSASpec(FreqVec,Para)

s11 = Para(1);
s22 = Para(2);
s33 = Para(3);
GaussianLB = Para(4);
LorentzianLB = Para(5);
Intensity = Para(6);

RefAxis = FreqVec;
AxisStart = 2*RefAxis(1) - RefAxis(2);
AxisIncrement = RefAxis(2) - RefAxis(1);
AxisLength = length(RefAxis);
FFTLength = AxisLength;

Spec  = zeros(1, FFTLength);
Spec1 = zeros(1, FFTLength);
Spec2 = zeros(1, FFTLength);

%
%define the apodization functions
%
  Qlb = abs(LorentzianLB);
  Qgb = abs(GaussianLB);
  QGaussB =  (Qgb/(sqrt(8*log(2))));
  Qmiddle = floor(FFTLength/2) + 1;
  Qemacht = fftshift(exp( -Qlb*(2*pi/(2*abs(AxisIncrement)*(FFTLength-1)))*abs(Qmiddle-(1:FFTLength)) - 2*((QGaussB)*(2*pi/(2*( abs(AxisIncrement)*(FFTLength-1) )))*abs(Qmiddle-(1:FFTLength))).^2 ));
  Qemacht = Qemacht / max(Qemacht);

%
%Now calculate the separate CSA tensor lineshapes
%

%
%First sort the principal values according to s33 > s22 > s11
%
  SortVec = sort([s11 s22 s33]);
  s11 = SortVec(1);
  s22 = SortVec(2);
  s33 = SortVec(3);


%
%Then determine the 4 areas that need to be distinguished to calculate the powder pattern
%

  %
  %This part is for w < s11
  %
  FreqVec1 = FreqVec(find(FreqVec < s11));
  PlaceVec1 = zeros(1, length(FreqVec1));


  %
  %This part is for s22 > w >= s11
  %
  FreqVec2 = FreqVec(find((FreqVec < s22) & (FreqVec >= s11)));
  PlaceVec2 = ellipke( ((FreqVec2-s11).*(s33-s22)) ./ ((s33-FreqVec2).*(s22-s11))  ) ./ ( sqrt(s33-FreqVec2) .* sqrt(s22-s11) .* pi );


  %
  %This part is for s22 == FreqVec (just in case the assymptote lies on a
  %grid point
  %
  FreqVec3 = FreqVec(find(FreqVec == s22));
  if ((s33 ~= s22) & (s22 ~= s11))
    PlaceVec3 = (ellipke( ((FreqVec3-0.001-s11).*(s33-s22)) ./ ((s33-FreqVec3+0.001).*(s22-s11))  ) ./ ( sqrt(s33-FreqVec3+0.001) .* sqrt(s22-s11) .* pi ) + ...
                 ellipke( ((s22-s11).*(s33-FreqVec3-0.001)) ./ ((s33-s22).*(FreqVec3+0.001-s11))  ) ./ ( sqrt(FreqVec3+0.001-s11) .* sqrt(s33-s22) .* pi ))/2;
  elseif (s33 ~= s22)
    PlaceVec3 = ellipke( ((s22-s11).*(s33-FreqVec3-0.001)) ./ ((s33-s22).*(FreqVec3+0.001-s11))  ) ./ ( sqrt(FreqVec3+0.001-s11) .* sqrt(s33-s22) .* pi );
  elseif (s22 ~= s11)
    PlaceVec3 = ellipke( ((FreqVec3-0.001-s11).*(s33-s22)) ./ ((s33-FreqVec3+0.001).*(s22-s11))  ) ./ ( sqrt(s33-FreqVec3+0.001) .* sqrt(s22-s11) .* pi );
  else
    PlaceVec3 = 1;
  end


  %
  %This part is for s33 >= w > s22
  %
  FreqVec4 = FreqVec(find((FreqVec <= s33) & (FreqVec > s22)));
  PlaceVec4 = ellipke( ((s22-s11).*(s33-FreqVec4)) ./ ((s33-s22).*(FreqVec4-s11))  ) ./ ( sqrt(FreqVec4-s11) .* sqrt(s33-s22) .* pi );


  %
  %This part is for w > s33
  %
  FreqVec5 = FreqVec(find(FreqVec > s33));
  PlaceVec5 = zeros(1, length(FreqVec5));


%
%Combine separate parts and integrate powder pattern always to the same
%integral per frequency resolution
%
  ret = [PlaceVec1 PlaceVec2 PlaceVec3 PlaceVec4 PlaceVec5];
  ret = ret / sum(ret) / abs(FreqVec(2)-FreqVec(1));

  Spec2 = ret;

  Spec2 = Spec2 * Intensity;
  Spec1 = ifft(fftshift(Spec2), FFTLength);
  Spec1 = Spec1 .* Qemacht;
  Spec1 = fftshift(fft(Spec1, FFTLength));

  Spec = Spec1;
 