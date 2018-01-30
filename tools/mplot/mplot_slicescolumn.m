function axh_v = mplot_slicescolumn(im3d,position,clim)
% function axh_v = mplot_slicescolumn(im3d,position,clim)
%

if isstruct(im3d)
    sz = size(im3d.r);
    if numel(sz) == 2; sz = [sz 1]; end

    axh_v = [];
    for k = 1:sz(3)
        im2d = zeros(sz(2),sz(1),3);
        im2d(:,:,1) = im3d.r(:,:,k)';
        im2d(:,:,2) = im3d.g(:,:,k)';
        im2d(:,:,3) = im3d.b(:,:,k)';
        im2d = im2d./repmat(max(im2d,[],3),[1 1 3]).*repmat(im3d.bright(:,:,k)',[1 1 3]);
        bottom = position.dbottom*(k-1);
        axh = axes('position',[position.left bottom position.width position.height]);
        axh_v = [axh_v; axh];
        imagesc(im2d)
    end
    set(axh_v,'YDir','normal')
    axis(axh_v,'equal','tight','off')
    set(axh_v,'CLim',clim)   
else
    sz = size(im3d);
    if numel(sz) == 2; sz = [sz 1]; end

    axh_v = [];
    for k = 1:sz(3)
        im2d = im3d(:,:,k);
        bottom = position.dbottom*(k-1);
        axh = axes('position',[position.left bottom position.width position.height]);
        axh_v = [axh_v; axh];
        imagesc(im2d')
    end
    set(axh_v,'YDir','normal')
    axis(axh_v,'equal','tight','off')
    colormap('gray')
    set(axh_v,'CLim',clim)
end

