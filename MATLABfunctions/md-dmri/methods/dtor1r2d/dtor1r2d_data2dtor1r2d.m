function dtor1r2d = dtor1r2d_data2dtor1r2d(s, xps, dtor1r2d_nodes)

[dtor1r2d_nx6,r1,r2] = dtor1r2d_nodes2nx6r1r2(dtor1r2d_nodes,xps);

k = exp(-xps.btomega*dtor1r2d_nx6').*(1-exp(-xps.tr*r1')).*exp(-xps.te*r2');
snorm = max(s);
w = snorm*lsqnonneg(k,s/snorm);

dtor1r2d = dtor1r2d_nodesw2dist(dtor1r2d_nodes,w);
dtor1r2d = dtor1r2d_sort(dtor1r2d);




