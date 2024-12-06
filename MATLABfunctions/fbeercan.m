function E = fbeercan(q,Delta,teta,R,L,D)

E = zeros(size(q));

%calculations

%Roots of Bessel equation J'm=0. Each row in the matrix are the
%10 first roots for a given m value, alfa10=0.
%There are 11 different m-values listed. 

alfa = [[ 0 3.8317   7.01559   10.1735   13.3237   16.4706 19.6159   22.7601   25.9037   29.0468   ];
[1.84118   5.33144   8.53632   11.706   14.8636  18.0155   21.1644   24.3113   27.4571   30.6019];  % J'1=0
[3.05424   6.70613   9.96947   13.1704   16.3475 19.5129   22.6716   25.826    28.9777   32.1273];  % J'2=0
[4.20119   8.01524   11.3459   14.5858   17.7887 20.9725   24.1449   27.3101   30.4703   33.6269];  % J'3=0
[5.31755   9.2824    12.6819   15.9641   19.196  19.196    22.401    25.5898   28.7678   31.9385];  % J'4=0
[6.41562  10.5199    13.9872   17.3128   20.5755 23.8036   27.0103   30.2028   33.3854   36.5608];  % J'5=0
[7.50127  11.7349    15.2682   18.6374   21.9317 25.1839   28.4098   31.6179   34.8134   37.9996];  % J'6=0
[8.57784  12.9324    16.5294   19.9419   23.2681 26.545    29.7907   33.0152   36.2244   39.4223];  % J'7=0
[9.64742  14.1155    17.774    21.2291   24.5872 27.8893   31.1553   34.3966   37.6201   40.8302];  % J'8=0
[10.7114  15.2867    19.0046   22.5014   25.8913 29.2186   32.5052   35.7638   39.0019   42.2246];  % J'9=0
[11.7709  16.4479    20.223    23.7607   27.182  30.5345   33.842    37.118    40.3711   43.6068]]; %J'10=0

nn = 5000;

for n = 1:nn
    alfakm2(:,:,n) = alfa.^2;
end

Size = size(alfa); mm = Size(1,1); kk = Size(1,2);

%Generate Knm matrix with Knm = 4 for n and m­0, 2 for one of n or m ­0 and 1 for 
%both n and m =0.

Knmmat = 4*ones(mm,kk,nn);
Knmmat(1,:,1)=1;
Knmmat(2:mm,:,1) = 2*ones(mm-1,kk,1);
Knmmat(1,:,2:nn) = 2*ones(1,kk,nn-1);



beta = 2*pi*q*R;

C1 = 2*R^2*beta.^4.*(sin(2*teta)).^2/L^2;
C2 = beta*L.*cos(teta)/R;
C3 = (pi*R/L)^2;
betac2 = (beta.*cos(teta)).^2;
betas = (beta.*sin(teta)); betas2 = (beta.*sin(teta)).^2;

m = 1:mm;
k = 1:kk;
n = 1:nn;

[m,k,n] = ndgrid(m,k,n);
onesmkn = ones(mm,kk,nn);


%Start by chosing a value of  beta
nmax = length(beta);
for b=1:nmax
    b
    x = betas(b)*onesmkn;
    x2 = betas2(b)*onesmkn;
    y = C2(b)*onesmkn;
    y2 = betac2(b)*onesmkn;

    %First, compute J'm(beta*sin(teta))
    Jprimm2 = (besselj(m-2,x) - (m-1)./x.*besselj(m-1,x)).^2;
    Jprimm2(1,:,:)=(besselj(1,x(1,1,1)))^2;


    %Now, compute "cylinder matrix"
    up=alfakm2.*Jprimm2.*exp(-alfakm2*D*Delta/R^2);          
    down = (alfakm2-x2).^2.*(alfakm2-(m-1).^2) + eps;
    cylmat = up./down;
    
    %Since alfa10=0 the entry 1,1 in cylmat matrix is of type 0/0
    %This entry is caculated separetly.
    up1 = Jprimm2(1,1,1)*exp(-alfakm2(1,1,1)*D*Delta/R^2);
    down1 = (alfakm2(1,1,1)-x2(1,1,1))^2;
    cylmat(1,1,:)=up1/down1;

    %Now, generate plane part
    upp=(1-((-1).^(n-1)).*cos(y)).*exp(-D*((n-1)*pi).^2*Delta/L^2);               
    downp = (C3.*(n-1).^2-y2).^2;
    planemat = upp./downp;

    %Put everything together
    totmat = Knmmat.*cylmat.*planemat;
    Sum1 = sum(totmat,1);
    Sum2 = sum(Sum1,2);
    Sum3 = sum(Sum2,3);

    E(b)  = C1(b)*Sum3;
    %End calculation for different beta values
end