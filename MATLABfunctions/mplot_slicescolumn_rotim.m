function [imh_v,axh_v] = mplot_slicescolumn_rotim(fh,im3d,position,clim)
% function axh_v = mplot_slicescolumn(im3d,position,clim)
%

if isstruct(im3d)
    sz = size(im3d.r);
    if numel(sz) == 2; sz = [sz 1]; end

    imh_v = [];
    axh_v = [];
    for k = 1:sz(3)
        left = position.papersize(1)*rand(1,1);
        bottom = position.papersize(2)*rand(1,1) - position.height/2;
        rotangle = 360*rand(1,1);
        magn = (1+1i)*exp(1i*rotangle/360*2*pi); magn = max(abs([real(magn) imag(magn)]));
        positionvector = [left bottom magn*[position.width position.height]];
        im2d = zeros(sz(2),sz(1),3);
        im2d(:,:,1) = im3d.r(:,:,k)';
        im2d(:,:,2) = im3d.g(:,:,k)';
        im2d(:,:,3) = im3d.b(:,:,k)';
        im2d = im2d./repmat(max(im2d,[],3),[1 1 3]).*repmat(im3d.bright(:,:,k)',[1 1 3]);
        axh = axes(fh,'Units','centimeters','position',positionvector);
        axh_v = [axh_v; axh];
        imh = imagesc(im2d);
        imh_v = [imh_v; imh];
        view(axh,[rotangle 90])
    end
    set(axh_v,'YDir','normal')
    axis(axh_v,'tight','off')
    set(axh_v,'CLim',clim)   
else
    sz = size(im3d);
    if numel(sz) == 2; sz = [sz 1]; end

    imh_v = [];
    axh_v = [];
    for k = 1:sz(3)
        left = position.papersize(1)*rand(1,1);
        bottom = position.papersize(2)*rand(1,1) - position.height/2;
        rotangle = 360*rand(1,1);
        magn = (1+1i)*exp(1i*rotangle/360*2*pi); magn = max(abs([real(magn) imag(magn)]));
        positionvector = [left bottom magn*[position.width position.height]];
        im2d = im3d(:,:,k);
        axh = axes(fh,'Units','centimeters','position',positionvector);
        imh = imagesc(im2d');        
        colormap(axh,gray(64))
        imh_v = [imh_v; imh];
        axh_v = [axh_v; axh];
        view(axh,[rotangle 90])
    end
    set(axh_v,'YDir','normal')
    axis(axh_v,'tight','off')
    set(axh_v,'CLim',clim)
end

