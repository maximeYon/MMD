function dtd_covariance_plot(S, xps, h, h2)
% function dtd_covariance_plot(S, xps, h, h2)

if (nargin < 4), h2 = []; end

opt = dtd_covariance_opt();
[m,cond,n] = dtd_covariance_1d_data2fit(S, xps, opt);
S_fit = dtd_codivide_1d_fit2data(m, xps);


plot(h, S);
title(h, num2str(n));




if (1)
    
    C = m(8:end) * 1e18;
    
    [x,y,z] = sphere(80);
    
    c = zeros(size(x));
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            r = [x(i,j) y(i,j) z(i,j)];
            
            r2 = tm_1x3_to_1x6(1,0,r);
            
            r4 = tm_1x6_to_1x21(tm_1x3_to_1x6(1,0,r));
            r4 = r4 / sqrt(tm_inner(r4, r4));
            L1 = tm_inner(C, r4);
            
            r4 = tm_1x6_to_1x21(tm_1x3_to_1x6(0,1,r));
            r4 = r4 / sqrt(tm_inner(r4, r4));
            L2 = tm_inner(C, r4);
            

            L = L1;
            
            x(i,j) = x(i,j) * L;
            y(i,j) = y(i,j) * L;
            z(i,j) = z(i,j) * L;
            
            c(i,j) = (L1 - L2);
            
        end
    end
    
    min(c(:))
    max(c(:))
    
    surface(h2,z,x,y,c);
    
    caxis(h2, [-0.5 0.5] );
    colormap(h2, fliplr(jet(100)));
    axis(h2, 'off')
    axis(h2, 'tight');
    axis(h2, 'equal');
    axis(h2, 'vis3d');
    camlight left;
    shading(h2, 'interp');
%     rotate3d(h2, 'on');
    material([0.7 0.4 0.8 1]); 

end