function [dpar,dperp,theta,phi,r1,r2,w] = my_dtr1r2d_bsmerge_initializepars(size_bsn,Nbsn)

dpar = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
dperp = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
theta = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
phi = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
r1 = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
r2 = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
w = single(zeros(size_bsn(1,1),size_bsn(1,2),size_bsn(1,3),size_bsn(1,4)*Nbsn));
