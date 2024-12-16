function mplot_technicolor_slices(method, dps, slices_path, clim, opt)
% function mplot_technicolor_slices(method, dps, slices_path, clim, opt)

if nargin < 5
    opt = mdm_opt();
end
opt = mplot_opt(opt);

sz = ones(1,3);
sz_temp = size(dps.s0);
sz(1:numel(sz_temp)) = sz_temp;

if isfield(opt,'k_range')
    k_range = opt.k_range;
else
    k_range = 1:sz(3); 
end

s0_nk = dps.s0(:,:,k_range);
s_b2000_nk = dps.s_b2000(:,:,k_range);
% smax = max(dps.s0(:,:,nk),[],'all');
% smax = quantile(s0_nk(isfinite(s0_nk)),.99,'all');
smax = quantile(s0_nk(s0_nk>0),.999,'all');
clim.s_b2000 = clim.s0;
clim.s0 = smax*clim.s0;
% clim.s2000 = max(dps.s2000(:,:,nk),[],'all')*clim.s2000;
% clim.s2000 = quantile(s2000_nk(isfinite(s0_nk)),.99,'all')*clim.s2000;
clim.s2000 = quantile(s_b2000_nk(s0_nk>0),.999,'all')*clim.s_b2000;

%Plot parameter maps
figure(2), clf

if strcmp(method,'dtr2d')
    plotfields.gray = {'s0';'s_b2000';'mdiso';'msddelta';'mr2';'vdiso';'vsddelta';'vr2'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor2';'cvsddeltar2'};
    plotfields.bin = {'mdiso';'msddelta';'mr2'};
elseif strcmp(method,'dtr1d')
    plotfields.gray = {'s0';'s_b2000';'mdiso';'msddelta';'mr1';'vdiso';'vsddelta';'vr1'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvsddeltar1'};
    plotfields.bin = {'mdiso';'msddelta';'mr1'};
elseif strcmp(method,'dtr1r2d')
    plotfields.gray = {'s0';'mdiso';'msddelta';'mr1';'mr2';'vdiso';'vsddelta';'vr1';'vr2'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvdisor2';'cvsddeltar1';'cvsddeltar2';'cvr1r2'};
    plotfields.bin = {'s0';'mdiso';'msddelta';'mr1';'mr2'};
elseif strcmp(method,'dtd')
    plotfields.gray = {'s0';'s_b2000';'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.hotcold = {'cvdisosddelta'};
    plotfields.bin = {'mdiso';'msddelta'};
elseif strcmp(method,'dtod')
    plotfields.gray = {'s0';'s_b2000';'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.hotcold = {'cvdisosddelta';'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu'};
%     plotfields.gray = {'s0';'s_b2000';'mdiso';'msddelta';'vdiso';'vsddelta';'dmdisodnu'};
%     plotfields.hotcold = {'cvdisosddelta';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu'};
    plotfields.bin = {'mdiso';'msddelta'};
    plotfields.omega.gray = {'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.omega.hotcold = {'cvdisosddelta'};
elseif strcmp(method,'dtor1r2d')
    plotfields.gray = {'s0';'s_b2000';'mdiso';'msddelta';'mr1';'mr2';'vdiso';'vsddelta';'vr1';'vr2'};
    plotfields.hotcold = {'cvdisosddelta';'cvdisor1';'cvdisor2';'cvsddeltar1';'cvsddeltar2';'cvr1r2';'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu'};
    plotfields.bin = {'mdiso';'msddelta';'mr1';'mr2'};
    plotfields.omega.gray = {'mdiso';'msddelta';'vdiso';'vsddelta'};
    plotfields.omega.hotcold = {'cvdisosddelta'};
elseif strcmp(method,'dtd_covariance')
    plotfields.gray = {'s0';'MD';'FA';'C_MD';'C_mu';'MKi';'MKa';'MKt';'MK'};
    plotfields.hotcold = {};
    plotfields.bin = {};
end

pixaspect = dps.nii_h.pixdim(3)/dps.nii_h.pixdim(2);
imaspect = sz(2)/sz(1);

% Nslices = numel(nk);
% Nparams = numel(plotfields.gray) + numel(plotfields.hotcold) + numel(dps.bin)*numel(plotfields.bin) + 4;
Nparams = 1;
Nslices = 1;
papersize = 3*[Nparams/opt.mplot.zoomx Nslices*pixaspect*imaspect/opt.mplot.zoomy];
dleft = 0;
pformat = '-dpdf';

position.left = 0;
position.dbottom = 0;
position.height = 1.01; 
position.width = 1.01; 

mask = dps.s0 > clim.mask_threshold*smax;
msf_mkdir(slices_path);
for nk = k_range
    msf_mkdir(fullfile(slices_path,num2str(nk)));

    for c = 1:numel(plotfields.gray)
        figure(2), clf
        im3d = dps.(plotfields.gray{c});
        position.left = position.left + dleft;
        axh_v = mplot_slicescolumn(mask(:,:,nk).*im3d(:,:,nk),position,clim.(plotfields.gray{c}));
        for n = 1:numel(axh_v), colormap(axh_v(n),gray(64)), end
        if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
            set(axh_v,'YDir','reverse')
        elseif strcmp(mdm_nii_oricode(dps.nii_h),'RAS')
            set(axh_v,'XDir','reverse','YDir','normal')
        else
            set(axh_v,'YDir','normal')
        end
        axis(axh_v,'tight','off')
        set(axh_v,'XLim',sz(1)/2 + 1/opt.mplot.zoomx*sz(1)/2*[-1 1],'YLim',sz(2)/2 + 1/opt.mplot.zoomy*sz(2)/2*[-1 1])

        set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
        set(gcf, 'InvertHardcopy','off','Color',[0 0 0]); 
        fig_fn = fullfile(slices_path,num2str(nk),plotfields.gray{c});
        print(fig_fn,'-loose',pformat)
    end

    for c = 1:numel(plotfields.hotcold)
        figure(2), clf
        im3d = dps.( plotfields.hotcold{c});
        position.left = position.left + dleft;
        axh_v = mplot_slicescolumn(mask(:,:,nk).*im3d(:,:,nk),position,clim.(plotfields.hotcold{c}));
        for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(64)), end
        if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
            set(axh_v,'YDir','reverse')
        elseif strcmp(mdm_nii_oricode(dps.nii_h),'RAS')
            set(axh_v,'XDir','reverse','YDir','normal')
        else
            set(axh_v,'YDir','normal')
        end
        axis(axh_v,'tight','off')
        set(axh_v,'XLim',sz(1)/2 + 1/opt.mplot.zoomx*sz(1)/2*[-1 1],'YLim',sz(2)/2 + 1/opt.mplot.zoomy*sz(2)/2*[-1 1])
        
        set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
        set(gcf, 'InvertHardcopy','off','Color',[0 0 0]); 
        fig_fn = fullfile(slices_path,num2str(nk),plotfields.hotcold{c});
        print(fig_fn,'-loose',pformat)
    end

    for c = 1:numel(plotfields.bin)
        for nbin = 1:numel(dps.bin)
            figure(2), clf
            clear im3d
            cind = (dps.bin{nbin}.(plotfields.bin{c})(:,:,nk)-min(clim.(plotfields.bin{c})))...
                /(max(clim.(plotfields.bin{c}))-min(clim.(plotfields.bin{c})));
            im3d = dist_cind2rgb_jet(cind);
            im3d.bright = mask(:,:,nk).*dps.bin{nbin}.f(:,:,nk);
            position.left = position.left + dleft;
            axh_v = mplot_slicescolumn(im3d,position,[0 1]);
            if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
                set(axh_v,'YDir','reverse')
            elseif strcmp(mdm_nii_oricode(dps.nii_h),'RAS')
                set(axh_v,'XDir','reverse','YDir','normal')
            else
                set(axh_v,'YDir','normal')
            end
            axis(axh_v,'tight','off')
            set(axh_v,'XLim',sz(1)/2 + 1/opt.mplot.zoomx*sz(1)/2*[-1 1],'YLim',sz(2)/2 + 1/opt.mplot.zoomy*sz(2)/2*[-1 1])
            
            set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
            set(gcf, 'InvertHardcopy','off','Color',[0 0 0]); 
            fig_fn = fullfile(slices_path,num2str(nk),[plotfields.bin{c} ['_bin' num2str(nbin)]]);
            print(fig_fn,'-loose',pformat)
        end
    end

    if isfield(dps,'bin')
        for nbin = 1:numel(dps.bin)
            figure(2), clf
            clear im3d
            im3d.r = dps.bin{nbin}.mdxx(:,:,nk);
            im3d.g = dps.bin{nbin}.mdyy(:,:,nk);
            im3d.b = dps.bin{nbin}.mdzz(:,:,nk);
            im3d.bright = mask(:,:,nk).*dps.bin{nbin}.f(:,:,nk);
            position.left = position.left + dleft;
            axh_v = mplot_slicescolumn(im3d,position,[0 1]);
            if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
                set(axh_v,'YDir','reverse')
            elseif strcmp(mdm_nii_oricode(dps.nii_h),'RAS')
                set(axh_v,'XDir','reverse','YDir','normal')
            else
                set(axh_v,'YDir','normal')
            end
            axis(axh_v,'tight','off')
            set(axh_v,'XLim',sz(1)/2 + 1/opt.mplot.zoomx*sz(1)/2*[-1 1],'YLim',sz(2)/2 + 1/opt.mplot.zoomy*sz(2)/2*[-1 1])
            
            set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
            set(gcf, 'InvertHardcopy','off','Color',[0 0 0]); 
            fig_fn = fullfile(slices_path,num2str(nk),['mdii' ['_bin' num2str(nbin)]]);
            print(fig_fn,'-loose',pformat)
        end
        
        figure(2), clf
        clear im3d
        im3d.r = dps.bin{1}.f(:,:,nk);
        im3d.g = dps.bin{2}.f(:,:,nk);
        im3d.b = dps.bin{3}.f(:,:,nk);
%         im3d.bright = mask(:,:,nk).*(dps.bin{1}.f(:,:,nk) + dps.bin{2}.f(:,:,nk) + dps.bin{3}.f(:,:,nk));
        im3d.bright = mask(:,:,nk);
        climRGB = [0 1];
        position.left = position.left + dleft;
        axh_v = mplot_slicescolumn(im3d,position,climRGB);
        if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
            set(axh_v,'YDir','reverse')
        elseif strcmp(mdm_nii_oricode(dps.nii_h),'RAS')
            set(axh_v,'XDir','reverse','YDir','normal')
        else
            set(axh_v,'YDir','normal')
        end
        axis(axh_v,'tight','off')
        set(axh_v,'XLim',sz(1)/2 + 1/opt.mplot.zoomx*sz(1)/2*[-1 1],'YLim',sz(2)/2 + 1/opt.mplot.zoomy*sz(2)/2*[-1 1])
        
        set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize);
        set(gcf, 'InvertHardcopy','off','Color',[0 0 0]); 
        fig_fn = fullfile(slices_path,num2str(nk),'fractionsRGB');
        print(fig_fn,'-loose',pformat)
    end

    if isfield(dps,'omega')
        for nomega = 1:numel(dps.omega)
            for c = 1:numel(plotfields.omega.gray)
                figure(2), clf
                im3d = dps.omega{nomega}.(plotfields.omega.gray{c});
                position.left = position.left + dleft;
                axh_v = mplot_slicescolumn(mask(:,:,nk).*im3d(:,:,nk),position,clim.(plotfields.omega.gray{c}));
                for n = 1:numel(axh_v), colormap(axh_v(n),gray(64)), end
                if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
                    set(axh_v,'YDir','reverse')
                elseif strcmp(mdm_nii_oricode(dps.nii_h),'RAS')
                    set(axh_v,'XDir','reverse','YDir','normal')
                else
                    set(axh_v,'YDir','normal')
                end
                axis(axh_v,'tight','off')
                set(axh_v,'XLim',sz(1)/2 + 1/opt.mplot.zoomx*sz(1)/2*[-1 1],'YLim',sz(2)/2 + 1/opt.mplot.zoomy*sz(2)/2*[-1 1])

                set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
                set(gcf, 'InvertHardcopy','off','Color',[0 0 0]); 
                fig_fn = fullfile(slices_path,num2str(nk),[plotfields.omega.gray{c} ['_omega' num2str(opt.dtod.maps_omega(nomega)/2/pi,3) 'Hz']]);
                print(fig_fn,'-loose',pformat)
                
                for c = 1:numel(plotfields.omega.hotcold)
                    figure(2), clf
                    im3d = dps.( plotfields.omega.hotcold{c});
                    position.left = position.left + dleft;
                    axh_v = mplot_slicescolumn(mask(:,:,nk).*im3d(:,:,nk),position,clim.(plotfields.omega.hotcold{c}));
                    for n = 1:numel(axh_v), colormap(axh_v(n),mplot_cmaphotcold(64)), end
                    if strcmp(mdm_nii_oricode(dps.nii_h),'LPS')
                        set(axh_v,'YDir','reverse')
                    elseif strcmp(mdm_nii_oricode(dps.nii_h),'RAS')
                        set(axh_v,'XDir','reverse','YDir','normal')
                    else
                        set(axh_v,'YDir','normal')
                    end
                    axis(axh_v,'tight','off')
                    set(axh_v,'XLim',sz(1)/2 + 1/opt.mplot.zoomx*sz(1)/2*[-1 1],'YLim',sz(2)/2 + 1/opt.mplot.zoomy*sz(2)/2*[-1 1])

                    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
                    set(gcf, 'InvertHardcopy','off','Color',[0 0 0]); 
                    fig_fn = fullfile(slices_path,num2str(nk),[plotfields.omega.hotcold{c} ['_omega' num2str(opt.dtod.maps_omega(nomega)/2/pi,3) 'Hz']]);
                    print(fig_fn,'-loose',pformat)
                end
            end
        end
    end

%     if ~isempty(slices_path)
%         msf_mkdir(fileparts(slices_path));
%         print(slices_path,'-loose',pformat)
%     end
end



