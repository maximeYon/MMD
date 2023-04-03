%
%Here we calculate the frequencies encountered for all powder averaging
%points in the mesh and for all the time steps on a rotor period at once.
%So, alpha and beta and t can be arrayed for optimum speed.
%This is eq. 4.17 from Matthias' 2006 lecture notes.
%
function Freq = FrequencyMAS(omega0, deltaaniso, eta, omegar, alpha, beta, gamma, t)

  Freq = zeros(length(alpha), length(t));

  alpha = alpha(:);
  beta = beta(:);
  gamma=gamma(:);
  t = t(:).';

  sinBeta = sin(beta);
  cosBeta = cos(beta);
  sin2Alpha = sin(2*alpha);
  cos2Alpha = cos(2*alpha);

  C1 = sqrt(2)/3 * omega0 * deltaaniso * sinBeta .* cosBeta .* (3 + eta*cos2Alpha);
  S1 = sqrt(2)/3 * omega0 * deltaaniso .* sinBeta * eta .* sin2Alpha;
  C2 = -omega0 * deltaaniso / 3 * ( 3/2*sinBeta.^2 - eta/2*(1+cosBeta.^2).*cos2Alpha );
  S2 =  omega0 * deltaaniso / 3 * eta * cosBeta .* sin2Alpha;

  if (length(gamma) == 1) 		%if gamma is a constant then use a fast vectorized approach
    Freq = C1 * cos(omegar*t - gamma) + S1 * sin(omegar*t-gamma) + C2 * cos(2*omegar*t - 2*gamma) + S2 * sin(2*omegar*t - 2*gamma);

  else 					%otherwise assume that the array in t is shorter than the array the angles
    for tel=1:length(t)
      Freq(:, tel) = C1 .* cos(omegar*t(tel) - gamma) + S1 .* sin(omegar*t(tel)-gamma) + C2 .* cos(2*omegar*t(tel) - 2*gamma) + S2 .* sin(2*omegar*t(tel) - 2*gamma);
    end
  end
