
function Pout = findminenergy(Pin,d)

global Niter energy

Niter = 1;

Pnorm = Pin;
%Pnorm = 1e-5;

options = optimset('Display', 'off','TolFun',1e-3,'TolX',1e-3,'MaxIter',1e9,'MaxFunEvals',1e9);
Pout = fminsearch(@fenergy, Pin./Pnorm, options);
Pout = Pout.*Pnorm;

 function energy = fenergy(Pin)
    Pin = Pin.*Pnorm;

    theta = Pin(1:d.N);
    phi = Pin((d.N+1):(2*d.N));

    x = sin(theta).*cos(phi);
    y = sin(theta).*sin(phi);
    z = cos(theta);

    [theta1,theta2] = ndgrid(theta,theta);
    [phi1,phi2] = ndgrid(phi,phi);
    x1 = sin(theta1).*cos(phi1);
    y1 = sin(theta1).*sin(phi1);
    z1 = cos(theta1);
    x2 = sin(theta2).*cos(phi2);
    y2 = sin(theta2).*sin(phi2);
    z2 = cos(theta2);

    r = sqrt(((x2-x1).^2+(y2-y1).^2+(z2-z1).^2));
    energymat = 1./r;
    diagelements = 1:(d.N+1):d.N^2;
    energymat(diagelements) = 0;
    energy = sum(reshape(energymat,numel(energymat),1));

    
    if any(Niter == 1e2:1e2:1e6)
        [Niter energy]
        figure(1), clf
        plot3(x,y,z,'o')
        axis(1.2*[-1 1 -1 1 -1 1])
        axis equal
        pause(.1)
    end
    Niter = Niter+1;
 end
end

