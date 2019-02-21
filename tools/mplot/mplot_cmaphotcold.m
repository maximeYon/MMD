function cmap = mplot_cmaphotcold(n)

cmap_hot = colormap(hot(n));
cmap_cold = [cmap_hot(:,3) cmap_hot(:,2) cmap_hot(:,1)];
cmap_hotcold = [flipud(cmap_cold); zeros(1,3); cmap_hot];

cmap = cmap_hotcold;
