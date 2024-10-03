clear all

M = 1000;
N = 120;

A = zeros(M,N);

tic

parfor n = 1:N
    E = eig(rand(M,M));
    A(:,n) = E;
end

toc
