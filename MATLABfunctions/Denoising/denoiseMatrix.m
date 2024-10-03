function [newX,sigma2,p] = denoiseMatrix(X)
% Takes as input matrix X with dimension MxN with N corresponding to the
% number of pixels and M to the number of data points. The output consists
% of "newX" containing a denoised version of X, "sigma2" an
% approximation to the data variation, "p" the number of signal carrying
% components.
[M,N] = size(X);
minMN = min([M N]);
% [U,S,V] = svd(X/sqrt(N),'econ');
[U,S,V] = svdecon(X/sqrt(N));
lambda = zeros(M,1);
lambda(1:minMN) = diag(S).^2;

p = 0;
pTest = false;
scaling = (M-(0:minMN))/N;
scaling(scaling<1) = 1;
while ~pTest
    sigma2 = (lambda(p+1)-lambda(minMN))/(4*sqrt((M-p)/N));
    pTest = sum(lambda(p+1:minMN))/scaling(p+1) >= (minMN-p)*sigma2;
    if ~pTest, p = p+1; end
end
sigma2 = sum(lambda(p+1:minMN))/(minMN-p)/scaling(p+1);

% newS = S;
% newS(p+1:end,p+1:end) = 0;
% newX = sqrt(N)*U*newS*V';
newX = sqrt(N)*U(:,1:p)*S(1:p,1:p)*V(:,1:p)';
end

function [U,S,V] = svdecon(X)
    % downloaded from: http://www.mathworks.com/matlabcentral/fileexchange/47132-fast-svd-and-pca/content/svdecon.m
    % Vipin Vijayan (2014)
    [m,n] = size(X);
    if  m <= n
        C = X*X';
        [U,D] = eig(C);
        clear C;

        [d,ix] = sort(abs(diag(D)),'descend');
        U = U(:,ix);    

        V = X'*U;
        s = sqrt(d);
        V = bsxfun(@(x,c)x./c, V, s');
        S = diag(s);
    else
        C = X'*X; 
        [V,D] = eig(C);
        clear C;

        [d,ix] = sort(abs(diag(D)),'descend');
        V = V(:,ix);    

        U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
        %s = sqrt(sum(U.^2,1))';
        s = sqrt(d);
        U = bsxfun(@(x,c)x./c, U, s');
        S = diag(s);
    end
end