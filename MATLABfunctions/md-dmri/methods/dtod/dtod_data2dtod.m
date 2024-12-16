function dtod = dtod_data2dtod(s, xps, dtod_nodes)

dtod_nx6 = dtod_nodes2nx6(dtod_nodes,xps);

k = exp(-xps.btomega*dtod_nx6');
snorm = max(s);
w = snorm*lsqnonneg(k,s/snorm);

dtod = dtod_nodesw2dist(dtod_nodes,w);
dtod = dtod_sort(dtod);




