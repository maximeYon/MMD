function mplot_technicolor_slices(method, dps, slices_path, clim, opt)
% function mplot_technicolor_slices(method, dps, slices_path, clim, opt)

%Plot parameter maps
figure(1), clf

smax = max(dps.s0(:));
clim.s0 = smax*clim.s0;

if strcmp(method,'dtr2d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr2';'vdiso';'vsddelta';'vr2'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor2';'cvsddeltar2'};
    plotfields.bin = {'mdiso';'msddelta';'mr2'};
elseif strcmp(method,'dtr1d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr1';'vdiso';'vsddelta';'vr1'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvsddeltar1'};
    plotfields.bin = {'mdiso';'msddelta';'mr1'};
elseif strcmp(method,'dtd')
    plotfields.gray = {'s0';'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.hotcold = {'cvdisosddelta'};
    plotfields.bin = {'mdiso';'msddelta'};
elseif strcmp(method,'dtd_covariance')
    plotfields.gray = {'s0';'MD';'FA';'C_MD';'C_mu';'MKi';'MKa';'MKt';'MK'};
    plotfields.hotcold = {};
    plotfields.bin = {};
end

sz = ones(1,3);
sz_temp = size(dps.s0);
sz(1:numel(sz_temp)) = sz_temp;
pixaspect = dps.nii_h.pixdim(3)/dps.nii_h.pixdim(2);
imaspect = sz(2)/sz(1);

% Nslices = numel(nk);
% Nparams = numel(plotfields.gray) + numel(plotfields.hotcold) + numel(dps.bin)*numel(plotfields.bin) + 4;
Nparams = 1;
Nslices = 1;
papersize = 3*[Nparams Nslices*pixaspect*imaspect];
dleft = 0;

position.left = 0;
position.dbottom = 0;
position.height = 1.01; 
position.width = 1.01; 

mask = dps.s0 > clim.mask_threshold*smax;
msf_mkdir(slices_path);
for nk = 1:sz(3) %nk = [2 6 11];
    msf_mkdir(fullfile(slices_path,num2str(nk)));

    for c = 1:numel(plotfields.gray)
        figure(1), clf
        im3d = dps.(plotfields.gray{c});
        position.left = position.left + dleft;
        axh_v = mplot_slicescolumn(mask(:,:,nk).*im3d(:,:,nk),position,clim.(plotfields.gray{c}));
        for n = 1:numel(axh_v), colormap(axh_v(n),gray(64)), end
        if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
            set(axh_v,'YDir','reverse')
        else
            set(axh_v,'YDir','normal')
        end

        set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
        fig_fn = fullfile(slices_path,num2str(nk),plotfields.gray{c});
        print(fig_fn,'-loose','-dpdf')
    end

    for c = 1:numel(plotfields.hotcold)
        figure(1), clf
        im3d = dps.( plotfields.hotcold{c});
        position.left = position.left + dleft;
        axh_v = mplot_slicescolumn(mask(:,:,nk).*im3d(:,:,nk),position,clim.(plotfields.hotcold{c}));
        for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(64)), end
        if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
            set(axh_v,'YDir','reverse')
        else
            set(axh_v,'YDir','normal')
        end
        
        set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
        fig_fn = fullfile(slices_path,num2str(nk),plotfields.hotcold{c});
        print(fig_fn,'-loose','-dpdf')
    end

    for c = 1:numel(plotfields.bin)
        for nbin = 1:numel(dps.bin)
            figure(1), clf
            clear im3d
            cind = (dps.bin{nbin}.(plotfields.bin{c})(:,:,nk)-min(clim.(plotfields.bin{c})))...
                /(max(clim.(plotfields.bin{c}))-min(clim.(plotfields.bin{c})));
            im3d = dist_cind2rgb_jet(cind);
            im3d.bright = mask(:,:,nk).*dps.bin{nbin}.f(:,:,nk);
            position.left = position.left + dleft;
            axh_v = mplot_slicescolumn(im3d,position,[0 1]);
            if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
                set(axh_v,'YDir','reverse')
            else
                set(axh_v,'YDir','normal')
            end
            
            set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
            fig_fn = fullfile(slices_path,num2str(nk),[plotfields.bin{c} num2str(nbin)]);
            print(fig_fn,'-loose','-dpdf')
        end
    end

    if isfield(dps,'bin')
        for nbin = 1:numel(dps.bin)
            figure(1), clf
            clear im3d
            im3d.r = dps.bin{nbin}.mdxx(:,:,nk);
            im3d.g = dps.bin{nbin}.mdyy(:,:,nk);
            im3d.b = dps.bin{nbin}.mdzz(:,:,nk);
            im3d.bright = mask(:,:,nk).*dps.bin{nbin}.f(:,:,nk);
            position.left = position.left + dleft;
            axh_v = mplot_slicescolumn(im3d,position,[0 1]);
            if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
                set(axh_v,'YDir','reverse')
            else
                set(axh_v,'YDir','normal')
            end
            
            set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
            fig_fn = fullfile(slices_path,num2str(nk),['mdii' num2str(nbin)]);
            print(fig_fn,'-loose','-dpdf')
        end
        
        figure(1), clf
        clear im3d
        im3d.r = dps.bin{1}.f(:,:,nk);
        im3d.g = dps.bin{2}.f(:,:,nk);
        im3d.b = dps.bin{3}.f(:,:,nk);
        im3d.bright = mask(:,:,nk).*(dps.bin{1}.f(:,:,nk) + dps.bin{2}.f(:,:,nk) + dps.bin{3}.f(:,:,nk));
        climRGB = [0 1];
        position.left = position.left + dleft;
        axh_v = mplot_slicescolumn(im3d,position,climRGB);
        if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
            set(axh_v,'YDir','reverse')
        else
            set(axh_v,'YDir','normal')
        end
        
        set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize);
        fig_fn = fullfile(slices_path,num2str(nk),'fractionsRGB');
        print(fig_fn,'-loose','-dpdf')
    end



%     if ~isempty(slices_path)
%         msf_mkdir(fileparts(slices_path));
%         print(slices_path,'-loose','-dpdf')
%     end
end



