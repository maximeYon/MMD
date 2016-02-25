function t = t2tpars(t3x3)
% function t = t2tpars(t3x3)
%
% calculate derived tensor parameters

t.t1x6 = tm_3x3_to_1x6(t3x3);
[V,D] = eig(t3x3);
lambdas = diag(D);
t.trace = sum(lambdas);
t.iso = t.trace/3;
[dummy,ind] = sort(lambdas,'descend');
t.lambdamax = lambdas(ind(1));                        
t.lambdamid = lambdas(ind(2));                       
t.lambdamin = lambdas(ind(3));  
t.lambdamaxvec = V(:,ind(1))';
t.lambdamidvec = V(:,ind(2))';
t.lambdaminvec = V(:,ind(3))';
Dlambdas = abs(lambdas-t.iso);
[dummy,ind] = sort(Dlambdas,'descend');
t.lambdazz = lambdas(ind(1));
t.lambdaxx = lambdas(ind(2));                        
t.lambdayy = lambdas(ind(3));
t.lambdazzvec = V(:,ind(1))';
t.lambdaxxvec = V(:,ind(2))';
t.lambdayyvec = V(:,ind(3))';

t.vlambda = 1/3*((t.lambdamax - t.iso).^2 + (t.lambdamid - t.iso).^2 + (t.lambdamin - t.iso).^2);
t.delta = (t.lambdazz - (t.lambdayy+t.lambdaxx)/2)./t.trace;
t.eta = 3*(t.lambdayy - t.lambdaxx+eps)./(2*t.trace.*t.delta+eps);
t.s = 3*t.lambdamin;
t.p = 2*(t.lambdamid - t.lambdamin);
t.l = t.lambdamax - t.lambdamid;
t.fa = sqrt(1/2)*sqrt((t.lambdamax-t.lambdamid).^2+(t.lambdamax-t.lambdamin).^2+(t.lambdamid-t.lambdamin).^2)...
                    ./sqrt(t.lambdamax.^2+t.lambdamid.^2+t.lambdamin.^2);
t.fa(isnan(t.fa)==1) = 0;
t.cl = (t.lambdamax - t.lambdamid)./t.lambdamax;
t.cl(isnan(t.cl)==1) = 0;
t.cp = (t.lambdamid - t.lambdamin)./t.lambdamax;
t.cp(isnan(t.cp)==1) = 0;

