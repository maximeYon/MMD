clear all

global Niter energy

for N = 51:200

    theta = acos(2*rand(N,1)-1)+2*pi;
    phi = 2*pi*rand(N,1)+2*pi;

    d.N = N;

    Pin = [theta; phi];
    Pout = findminenergy(Pin,d); Pin = Pout;
    theta = Pin(1:N);
    phi = Pin((N+1):(2*N));

    theta = angle(exp(i*theta));
    phi = angle(exp(i*phi));

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
    diagelements = 1:(N+1):N^2;
    energymat(diagelements) = 0;
    energy = sum(reshape(energymat,numel(energymat),1))

    figure(1), clf
    plot3(x,y,z,'o')
    axis(1.2*[-1 1 -1 1 -1 1])
    axis equal
    pause(1)

    eval(['save UniformDistSphereRepulsionN' num2str(N) ' theta phi'])

end