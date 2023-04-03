function s = fexsyr1r2_fit2data(m, xps)
% function s = fexsyr1r2_fit2data(m, xps)


% convert to readable parameters
meqi          = m(1);
meqe          = m(2);
R1i           = m(3);
R1e           = m(4);
R2i           = m(5);
R2e           = m(6);
Di            = m(7);
De            = m(8);
kie           = m(9);


kei = kie*meqi/meqe;

Meq = [meqi; meqe];
K = [ kie -kei
     -kie  kei];
R1 = diag([R1i; R1e]);
R2 = diag([R2i; R2e]);
D  = diag([ Di;  De]);


s = zeros(xps.n,1);

for ntd1 = 1:xps.n

    vtau_KR1 = xps.vtau_KR1(ntd1);
    vtau_KR2 = xps.vtau_KR2(ntd1);
    q1 = xps.q1(ntd1);
    q2 = xps.q2(ntd1);
    tau_KR2 = xps.tau_KR2(ntd1);
    tau_KR2D = xps.tau_KR2D(ntd1);
    tau_KR1D = xps.tau_KR1D(ntd1);
    tr = xps.tr(ntd1);

    blocks_Meq= [1 0     0 0 0  0 0 0     0 0 0  0 0]';
    blocks_K  = [1 1     1 1 1  1 1 1     1 1 1  1 1]';
    blocks_R1 = [1 0     0 1 0  0 1 0     0 1 0  0 0]';
    blocks_R2 = [0 1     1 0 1  1 0 1     1 0 1  1 1]';
    blocks_q  = [0 0 q1*[1 1 1] 0 0 0 q2*[1 1 1] 0 0]';
    blocks_tau = [tr       tau_KR2  ...
                  tau_KR2D tau_KR1D tau_KR2D ...
                  tau_KR2  vtau_KR1 tau_KR2 ...
                  tau_KR2D tau_KR1D tau_KR2D ...
                  tau_KR2  vtau_KR2]';
    blocks_n = numel(blocks_R1);

    M0 = [0; 0];
    M = M0;
    M_array = repmat(M,[1 blocks_n+1]);

    for nblock = 1:blocks_n
        A = -(blocks_K(nblock)*K + blocks_R1(nblock)*R1 + blocks_R2(nblock)*R2 + blocks_q(nblock)^2*D);
        [A_eigenvecs,A_eigenvals] = eig(A);
        propagator = A_eigenvecs*diag(exp(diag(A_eigenvals)*blocks_tau(nblock)))*inv(A_eigenvecs);
        M = propagator*(M - blocks_Meq(nblock)*Meq) + blocks_Meq(nblock)*Meq;
        M_array(:,nblock+1) = M;    
    end
% M_array
%     figure(1), clf
%     semilogy(cumsum([0; blocks_tau],1),M_array','-')
%     %set(gca,'YLim',sum(Meq)*[-.1 1.1])
%     title(num2str(ntd1))
%     pause(1)

    s(ntd1) = sum(M);
end
