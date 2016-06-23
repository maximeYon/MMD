function mplot_tensors_in_voxel(dt, n_tensor, sc)
% function mplot_tensors_in_voxel(dt, n_tensor, sc)
%
% WIP
%
% n_tensor - number of tensors to display per voxel side
% sc       - scaling factor (adjust to rescale tensors)

if (nargin < 2), n_tensor = 4; end
if (nargin < 3), sc = 0.1; end

do_plot_box = 1;
do_show_tensors = 1;



% plot box
if (do_plot_box)
    
    lcol = [0 0 0] + 0.8;
    a = 0.65;
    b = n_tensor + 1 - a;
    
    plot3([a b b a a], [a a b b a], [a a a a a ], 'k', 'color', lcol); hold on;
    plot3([a b b a a], [a a b b a], [b b b b b ], 'k', 'color', lcol);
    plot3([a a], [a a], [a b], 'k', 'color', lcol);
    plot3([a a], [b b], [a b], 'k', 'color', lcol);
    plot3([b b], [b b], [a b], 'k', 'color', lcol);
    plot3([b b], [a a], [a b], 'k', 'color', lcol);
    
end

if (do_show_tensors)
    [x,y,z] = sphere(29);
    
    
    for i = 1:n_tensor
        for j = 1:n_tensor
            for k = 1:n_tensor
              
                % select a tensor on random
                ind = randi(size(dt,1));
                
                dt_3x3 = tm_1x6_to_3x3(dt(ind,:));
                
                [M, eig_vals] = eigs(dt_3x3);
                eig_vals = diag(eig_vals) * 1e9;
                
                f = @(x,y,z) sc * [x(:) y(:) z(:)] * M;
                h = @(x) x(:);
                g = @(x,y,z) cellfun(@(p) reshape(p,size(x,1),size(x,2)), ...
                    mat2cell(f(x,y,z),numel(x), [1 1 1]), 'uniformoutput', 0);
                
                p = g(x * eig_vals(1),y * eig_vals(2),z * eig_vals(3));
                [x3,y3,z3] = deal(p{:});
                
                c = mean(eig_vals);
                
                r = randn(1,3) * 0.1;
                
                surf(x3 + i + r(1),y3 + j + r(2), z3 + k + r(3),zeros(size(x3)) + c); hold on;
            end
        end
    end
end

caxis([-1 3]);
colormap(fliplr(jet(100)));
axis off tight vis3d;
camlight left;
shading interp;
material([0.7 0.4 0.8 1] );