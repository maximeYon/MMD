function dtr2d = dtr2d_data2dtr2d(stemp, bt_mx6, te, dtr2d_nodes,opt)

r2reg = opt.dtr2d.r2reg;

[dtr2d_nx6,r2] = dtr2d_nodes2nx6r2(dtr2d_nodes);

k = exp(-bt_mx6*dtr2d_nx6').*exp(-te*r2');
maxk = max(k,[],1);
A = [k; r2reg*numel(stemp)*(1-maxk)];
snorm = max(stemp);
y = [stemp/snorm; 0];
    
w = snorm*lsqnonneg(A,y);

% ycalc = A*(w/snorm);
% chisq_v = (y - ycalc).^2;
% chisq_fit = sum(chisq_v(1:(end-1)));
% chisq_reg = chisq_v(end);
% [chisq_fit chisq_reg]

dtr2d = dtr2d_nodesw2dist(dtr2d_nodes,w);
dtr2d = dtr2d_sort(dtr2d);




