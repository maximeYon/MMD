function dtr2d = dtr2d_data2dtr2d(stemp, bt_mx6, te, dtr2d_nodes,opt)

r2reg = 0;

[dtr2d_nx6,r2] = dtr2d_nodes2nx6r2(dtr2d_nodes);

k = exp(-bt_mx6*dtr2d_nx6').*exp(-te*r2');
maxk = max(k,[],1);
A = [k; r2reg*(1-maxk)];
snorm = max(stemp);
y = [stemp/snorm; 0];
    
w = snorm*lsqnonneg(A,y);

dtr2d = dtr2d_nodesw2dist(dtr2d_nodes,w);
dtr2d = dtr2d_sort(dtr2d);




