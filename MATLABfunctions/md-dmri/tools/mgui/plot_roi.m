function [Ic, h_plot] = plot_roi(I_roi, col, lw, scale, h_axes)
% function plot_roi(I_roi, col, lw)

if (nargin < 2), col = 'r'; end
if (nargin < 3), lw = 1; end
if (nargin < 4), scale = 1; end
if (nargin < 5), h_axes = gca; end
if (sum(I_roi(:)) == 0), return; end

% Detect vertical and horizontal edges
ti = [];
tj = [];
for c_edge = 1:2
    switch (c_edge)
        case 1
            shift = [0 1];
            di    = [0 1];
            dj    = [0 0];
        case 2
            shift = [1 0];
            di    = [0 0];
            dj    = [0 1];
    end

    % Find edges
    I_edge = I_roi - circshift(I_roi, shift);

    % Turn edge indices into coordinates
    ind = find(I_edge ~= 0);
    [i,j] = ind2sub(size(I_roi), ind);
    i = i - 0.5;
    j = j - 0.5;

    % Plot each detected edge
    tj = [tj; ([j j] + repmat(dj, size(j,1), 1)) * scale];
    ti = [ti; ([i i] + repmat(di, size(i,1), 1)) * scale];
end

if (1)
    line(tj', ti', 'Color', col', 'Parent', h_axes, 'LineWidth', lw);
else
    plot(h_axes, tj', ti', '-', 'Color', col', 'LineWidth', lw);
end


% warning('Try using line instead of plot!');