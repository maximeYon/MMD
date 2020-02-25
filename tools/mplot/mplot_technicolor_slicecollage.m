function mplot_technicolor_slicecollage(method, dps, slices_path, clim, opt)
% function mplot_technicolor_slicecollage(method, dps, slices_path, clim, opt)

%Plot parameter maps

smax = max(dps.s0(:));
clim.s0 = smax*clim.s0;
clim.s2000 = max(dps.s2000(:))*clim.s2000;
mask = dps.s0 > clim.mask_threshold*smax;

sz = ones(1,3);
sz_temp = size(dps.s0);
sz(1:numel(sz_temp)) = sz_temp;
pixaspect = dps.nii_h.pixdim(3)/dps.nii_h.pixdim(2);
imaspect = sz(2)/sz(1);

if strcmp(method,'dtr2d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr2';'vdiso';'vsddelta';'vr2'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor2';'cvsddeltar2'};
    plotfields.bin = {'mdiso';'msddelta';'mr2'};
elseif strcmp(method,'dtr1d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr1';'vdiso';'vsddelta';'vr1'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvsddeltar1'};
    plotfields.bin = {'mdiso';'msddelta';'mr1'};
elseif strcmp(method,'dtd')
    plotfields.gray = {'s0';'s2000';'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.hotcold = {'cvdisosddelta'};
    plotfields.bin = {'mdiso';'msddelta'};
    Nparams = numel(plotfields.gray) + numel(plotfields.hotcold) + numel(dps.bin)*numel(plotfields.bin) + 4;
    papersize = 3*[7 3*pixaspect*imaspect];

    dbottom = 1/3;
    dleft = 1/7;

    position.height = 1.01*dbottom; 
    position.width = 1.01*dleft;
    
    position.left_v =   [0 0 5 6 5 6 5.5 1.5 1.5 1.5 2.5 2.5 2.5 3.5 3.5 3.5 0]/7;
    position.bottom_v = [2 1 2 2 1 1 0 2 1 0 2 1 0 2 1 0 0]/3;
end


msf_mkdir(fullfile(slices_path));
%msf_mkdir(fullfile(slices_path,'matfigs'));
for nk = 1:sz(3)
% for nk = 2
    figure(1), clf
    map_count = 1;
    axh_v = [];
    for c = 1:numel(plotfields.gray)
        im3d = double(dps.(plotfields.gray{c}));
        im2d = mask(:,:,nk).*im3d(:,:,nk);
        position.left = position.left_v(map_count);
        position.bottom = position.bottom_v(map_count);
        axh = axes('position',[position.left position.bottom position.width position.height]);
        imagesc(im2d')
        colormap(axh,gray(256))
        set(axh,'CLim',clim.(plotfields.gray{c}))
        axh_v = [axh_v; axh];
        map_count = map_count+1;
    end

    for c = 1:numel(plotfields.hotcold)
        im3d = double(dps.(plotfields.hotcold{c}));
        im2d = mask(:,:,nk).*im3d(:,:,nk);
        position.left = position.left_v(map_count);
        position.bottom = position.bottom_v(map_count);
        axh = axes('position',[position.left position.bottom position.width position.height]);

        imagesc(im2d')
        colormap(axh,mplot_cmaphotcold(128))
        set(axh,'CLim',clim.(plotfields.hotcold{c}))   
        axh_v = [axh_v; axh];
        map_count = map_count+1;

    end

    
    for c = 1:numel(plotfields.bin)
        for nbin = 1:numel(dps.bin)
            clear im3d
            cind = (dps.bin{nbin}.(plotfields.bin{c})-min(clim.(plotfields.bin{c})))...
                /(max(clim.(plotfields.bin{c}))-min(clim.(plotfields.bin{c})));
            im3d = dist_cind2rgb_jet(cind);
            im3d.bright = mask.*dps.bin{nbin}.f;
            im2d = zeros(sz(2),sz(1),3);
            im2d(:,:,1) = im3d.r(:,:,nk)';
            im2d(:,:,2) = im3d.g(:,:,nk)';
            im2d(:,:,3) = im3d.b(:,:,nk)';
            im2d = im2d./repmat(max(im2d,[],3),[1 1 3]).*repmat(im3d.bright(:,:,nk)',[1 1 3]);
            position.left = position.left_v(map_count);
            position.bottom = position.bottom_v(map_count);
            axh = axes('position',[position.left position.bottom position.width position.height]);
            
            imagesc(im2d)
            set(axh,'CLim',[0 1])   
            axh_v = [axh_v; axh];
            map_count = map_count+1;
        end
    end

    for nbin = 1:numel(dps.bin)
        clear im3d
        im3d.r = dps.bin{nbin}.mdxx;
        im3d.g = dps.bin{nbin}.mdyy;
        im3d.b = dps.bin{nbin}.mdzz;
        im3d.bright = mask.*dps.bin{nbin}.f;
        im2d = zeros(sz(2),sz(1),3);
        im2d(:,:,1) = im3d.r(:,:,nk)';
        im2d(:,:,2) = im3d.g(:,:,nk)';
        im2d(:,:,3) = im3d.b(:,:,nk)';
        im2d = im2d./repmat(max(im2d,[],3),[1 1 3]).*repmat(im3d.bright(:,:,nk)',[1 1 3]);
        position.left = position.left_v(map_count);
        position.bottom = position.bottom_v(map_count);
        axh = axes('position',[position.left position.bottom position.width position.height]);

        imagesc(im2d)
        set(axh,'CLim',[0 1])   
        axh_v = [axh_v; axh];
        map_count = map_count+1;
    end

    clear im3d
    im3d.r = dps.bin{1}.f;
    im3d.g = dps.bin{2}.f;
    im3d.b = dps.bin{3}.f;
    im3d.bright = mask.*(dps.bin{1}.f + dps.bin{2}.f + dps.bin{3}.f);
    im2d = zeros(sz(2),sz(1),3);
    im2d(:,:,1) = im3d.r(:,:,nk)';
    im2d(:,:,2) = im3d.g(:,:,nk)';
    im2d(:,:,3) = im3d.b(:,:,nk)';
    im2d = im2d./repmat(max(im2d,[],3),[1 1 3]).*repmat(im3d.bright(:,:,nk)',[1 1 3]);
    position.left = position.left_v(map_count);
    position.bottom = position.bottom_v(map_count);
    axh = axes('position',[position.left position.bottom position.width position.height]);

    imagesc(im2d)
    set(axh,'CLim',[0 1])   
    axh_v = [axh_v; axh];
    map_count = map_count+1;

    if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
        set(axh_v,'YDir','reverse')
    else
        set(axh_v,'YDir','normal')
    end
    axis(axh_v,'tight','off')

    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    fig_name = ['slice' num2str(nk)];
    print(fullfile(slices_path,fig_name),'-loose','-dpdf')
    %print(fullfile(slices_path,fig_name),'-loose','-dpng')
    %savefig(gcf,fullfile(slices_path,fig_name),'compact')
end
