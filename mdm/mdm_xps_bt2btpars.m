function xps_out = mdm_xps_bt2btpars(xps_in)
% function xps_out = mdm_xps_bt2btpars(xps_in)
%
% calculate derived b-tensor parameters

xps = xps_in;

% if isfield(xps,'b_delta')
%     error('xps field b_delta already exists')
% end
    
lambda.min = zeros(xps.n,1);
lambda.mid = zeros(xps.n,1);
lambda.max = zeros(xps.n,1);
xps.b_pasxx = zeros(xps.n,1); 
xps.b_pasyy = zeros(xps.n,1); 
xps.b_paszz = zeros(xps.n,1); 
xps.b_pasxxvec = zeros(xps.n,3); 
xps.b_pasyyvec = zeros(xps.n,3); 
xps.b_paszzvec = zeros(xps.n,3); 

for n = 1:xps.n
    bt = tm_1x6_to_3x3(xps.bt(n,:));
    [V,D] = eig(bt);
    lambdas = diag(D);
    [dummy,indx] = sort(lambdas,'descend');
    lambda.max(n,1) = min(lambdas(indx(1)));                        
    lambda.mid(n,1) = min(lambdas(indx(2)));                        
    lambda.min(n,1) = min(lambdas(indx(3)));    
    Dlambdas = abs(lambdas-sum(lambdas)/3);
    [dummy,indx] = sort(Dlambdas,'descend');
    xps.b_paszz(n,1) = lambdas(indx(1));
    xps.b_pasxx(n,1) = lambdas(indx(2));                        
    xps.b_pasyy(n,1) = lambdas(indx(3));
    xps.b_paszzvec(n,:) = V(:,indx(1))';
    xps.b_pasxxvec(n,:) = V(:,indx(2))';
    xps.b_pasyyvec(n,:) = V(:,indx(3))';
end

xps.b = xps.b_pasxx + xps.b_pasyy + xps.b_paszz;
xps.b_delta = (xps.b_paszz - (xps.b_pasyy+xps.b_pasxx)/2)./xps.b;
xps.b_eta = 3*(xps.b_pasyy - xps.b_pasxx+eps)./(2*xps.b.*xps.b_delta+eps);
xps.b_s = 3*lambda.min;
xps.b_p = 2*(lambda.mid - lambda.min);
xps.b_l = lambda.max - lambda.mid;

xps.u = xps.b_paszzvec;
xps.b_theta = acos(xps.u(:,3));
xps.b_phi = atan2(xps.u(:,2),xps.u(:,1));

%to be implemented
%xps.b_alpha   = ; 
%xps.b_beta   = ; 
%xps.b_gamma   = ; 

xps_out = xps;