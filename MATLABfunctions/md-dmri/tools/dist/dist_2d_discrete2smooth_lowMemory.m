function dist_s = dist_2d_discrete2smooth_lowMemory(dist_d,dist_s)
% function dist_s = dist_2d_discrete2smooth(dist_d,dist_s)
%
% Convert from discrete to smooth 2D distribution. Gaussian convolution.
%
% dist_d discrete distribution  with fields
% dist_d.n number of discrete components
% dist_d.x vector of x-values
% dist_d.y vector of y-values
% dist_d.w vector of weights
% dist_s smooth distribution  with input fields
% dist_s.x vector of x-values
% dist_s.y vector of y-values
% dist_s.xsigma std for Gaussian convolution x
% dist_s.ysigma std for Gaussian convolution y
% dist_s output smooth distribution with additional fields
% dist_s.w matrix of weights
% dist_s.wprojx vector of weights projected on x-axis
% dist_s.wprojy vector of weights projected on y-axis

if isstruct(dist_d)
    n0 = dist_d.n;
    x0 = dist_d.x;
    y0 = dist_d.y;
    w0 = dist_d.w;
else
    x0 = dist_d(:,1);
    y0 = dist_d(:,2);
    w0 = dist_d(:,3);
    n0 = numel(w0);
end

x = dist_s.x;
y = dist_s.y;
xsigma = dist_s.xsigma;
ysigma = dist_s.ysigma;
nx = numel(x);
ny = numel(y);

dx = x(2)-x(1);
dy = y(2)-y(1);
[xx,yy] = ndgrid(x,y);

%% low memory implementation
w0_2D = zeros(size(x,1),size(y,1));
% avoid folding of out of bound data
w0(x0>max(x))=0; w0(x0<min(x))=0;
w0(y0>max(y))=0; w0(y0<min(y))=0;
% loop
for ind = 1:n0
    [~,ind_x] = min(abs(x-x0(ind,1)));
    [~,ind_y] = min(abs(y-y0(ind,1)));
    w0_2D(ind_x,ind_y) = w0_2D(ind_x,ind_y) + w0(ind,1);
end
% my_Gaussian = fspecial('gaussian',nx,(xsigma./dx));
% dist_s.w = conv2(w0_2D,my_Gaussian,'same');
dist_s.w=imgaussfilt(w0_2D,[xsigma./dx ysigma./dy],'Padding',0);

% xx_k = repmat(xx(:),[1 n0]);
% yy_k = repmat(yy(:),[1 n0]);
% x0_k = repmat(x0',[nx*ny 1]);
% y0_k = repmat(y0',[nx*ny 1]);
% k = 1/(xsigma*sqrt(2*pi)).*exp(-(xx_k-x0_k).^2/(2*xsigma^2)).*...
%     1/(ysigma*sqrt(2*pi)).*exp(-(yy_k-y0_k).^2/(2*ysigma^2));
% 
% p = k*w0;
% dist_s.w = reshape(p,[nx ny]);

% figure(1)
% subplot(1,2,1)
% imagesc(dist_s.w)
% subplot(1,2,2)
% imagesc(w0_2D)

dist_s.wprojx = sum(dist_s.w,2);
dist_s.wprojy = sum(dist_s.w,1);