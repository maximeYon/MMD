function dtd = dtd_data2dtd(stemp,bt_mx6,dtd_nodes)

dtd_nx6 = dtd_nodes2nx6(dtd_nodes);

k = exp(-bt_mx6*dtd_nx6');
snorm = max(stemp);
w = snorm*lsqnonneg(k,stemp/snorm);

dtd = dtd_nodesw2dist(dtd_nodes,w);
dtd = dtd_sort(dtd);




