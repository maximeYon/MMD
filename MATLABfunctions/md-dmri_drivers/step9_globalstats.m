clear all

data_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419';


case_dir = dir(fullfile(data_path,'P*'));
case_names = cell(size(case_dir,1),1);
for ncase = 1:size(case_dir,1)
    case_names{ncase} = case_dir(ncase).name;
end

% case_names = {'P0275752'; 'P0246329'; 'P0254855'; 'P0262819'; 'P0286389'; 'P0286574'};
% case_names = {'P0066871'};

nii_ext = '.nii.gz';
nii_name = 'dtd_s2000';

for ncase = 1:numel(case_names)
    case_name = case_names{ncase};
    
    nii_path = fullfile(data_path,case_name,'rwi','pmaps');    
    bs_path = fullfile(fileparts(nii_path),'dtd','bootstraps');

    nii_fn = fullfile(nii_path,[nii_name nii_ext]);
    
    if exist(nii_path,'dir') ~= 7, disp([case_name ' nii_path missing']), continue, end
    if exist(bs_path,'dir') ~= 7, disp([case_name ' bs_path missing']), continue, end
    if exist(nii_fn,'file') ~= 2, disp([case_name ' nii_fn missing']), continue, end
        
    opt = mdm_opt();
    opt = mplot_opt(opt);
    opt.roiplot.case_name = case_name;

    [I,nii_h] = mdm_nii_read(nii_fn);

    %Plot distributions for ROI

    clim.mdiso = [0 4]*1e-9;
    clim.msddelta = [0 1];

    axpars.xmin = clim.mdiso(1);
    axpars.xmax = clim.mdiso(2);
    axpars.ymin = clim.msddelta(1);
    axpars.ymax = clim.msddelta(2);
    contourpars.Nx = 50;
    contourpars.Ny = contourpars.Nx;

    dist_s.x = linspace(axpars.xmin,axpars.xmax,contourpars.Nx)';
    dist_s.y = linspace(axpars.ymin,axpars.ymax,contourpars.Ny)';
    dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));
    dist_s.ysigma = 2*abs(dist_s.y(2) - dist_s.y(1));
    dist_s_in = dist_s;

    bsno = msf_getdirno(bs_path);        
    for nbs = 1:numel(bsno)
        mfs_fn   = fullfile(bs_path,num2str(bsno(nbs)),'mfs.mat');
        if exist(mfs_fn,'file')==2

            mfs = mdm_mfs_load(mfs_fn);
            sz = size(mfs.m);
            
            [dpar,dperp,~,~,w] = dtd_4d_m2pars(mfs.m);
            [diso,~,~,~,~,sddelta] = dtd_pars2dpars(dpar,dperp); 

            ind = logical(repmat(mfs.mask,[1 1 1 sz(4)]));
            break
        end
    end
    
    bs_w = NaN*ones([contourpars.Nx,contourpars.Ny,numel(bsno)]);
    tic
    parfor nbs = 1:numel(bsno)
        mfs_fn   = fullfile(bs_path,num2str(bsno(nbs)),'mfs.mat');
        if exist(mfs_fn,'file')==2
            bs_w_nk = 0*ones([contourpars.Nx,contourpars.Ny,sz(3)]);

            mfs = mdm_mfs_load(mfs_fn);
            [dpar,dperp,~,~,w] = dtd_4d_m2pars(mfs.m);
            [diso,~,~,~,~,sddelta] = dtd_pars2dpars(dpar,dperp);           
            
            for nk = 1:sz(3)
                diso_nk = diso(:,:,nk);
                sddelta_nk = sddelta(:,:,nk);
                w_nk = w(:,:,nk);
                ind_nk = ind(:,:,nk);
                dist_s_temp = dist_2d_discrete2smooth([diso_nk(ind_nk) sddelta_nk(ind_nk) w_nk(ind_nk)],dist_s_in);
                bs_w_nk(:,:,nk) = dist_s_temp.w;
            end
            bs_w(:,:,nbs) = sum(bs_w_nk,3);                
        end
    end
    toc

    dist_s.w = nanmedian(bs_w,3);
    diso_sddelta_dist_s = dist_s;
%     clear bs_w
    %%
    
    figure(2), clf
    
%     sz = size(I);
% 
%     clim.s2000 = 1*[0 1];
%     fov.x = nii_h.dim(2)*nii_h.pixdim(2);
%     fov.y = nii_h.dim(3)*nii_h.pixdim(3);
%     fov.z = nii_h.dim(4)*nii_h.pixdim(4);
% 
%     roi_projx = squeeze(sum(sum(roi_I,2),3));
%     roi_maxx = max(find(roi_projx==max(roi_projx)));
%     im2d_x = squeeze(I(roi_maxx,:,:));
%     roi2d_x = squeeze(roi_I(roi_maxx,:,:));
% 
%     roi_projy = squeeze(sum(sum(roi_I,1),3));
%     roi_maxy = max(find(roi_projy==max(roi_projy)));
%     im2d_y = squeeze(I(:,roi_maxy,:));
%     roi2d_y = squeeze(roi_I(:,roi_maxy,:));
% 
%     roi_projz = squeeze(sum(sum(roi_I,1),2));
%     roi_maxz = max(find(roi_projz==max(roi_projz)));
%     im2d_z = squeeze(I(:,:,roi_maxz));
%     roi2d_z = squeeze(roi_I(:,:,roi_maxz));
% 
%     axh_pmap_z = axes('position',[0 fov.z/(fov.y+fov.z) fov.x/(fov.x+fov.y) fov.y/(fov.y+fov.z)]);
%     imagesc(axh_pmap_z,im2d_z')
%     set(axh_pmap_z,'YDir','normal')
%     plot_roi(roi2d_z');
% 
%     axh_pmap_x = axes('position',[fov.x/(fov.x+fov.y) 0 fov.y/(fov.x+fov.y) fov.z/(fov.y+fov.z)]);
%     imagesc(axh_pmap_x,im2d_x')
%     set(axh_pmap_x,'YDir','normal')
%     plot_roi(roi2d_x');
% 
%     axh_pmap_y = axes('position',[0 0 fov.x/(fov.x+fov.y) fov.z/(fov.y+fov.z)]);
%     imagesc(axh_pmap_y,im2d_y')
%     set(axh_pmap_y,'YDir','normal')
%     plot_roi(roi2d_y');
% 
%     hold(axh_pmap_z,'on')
%     phxz = plot(axh_pmap_z,roi_maxx*[1 1],[-1 sz(1)+1],'-');
%     phyz = plot(axh_pmap_z,[-1 sz(2)+1],roi_maxy*[1 1],'-');
%     hold(axh_pmap_y,'on')
%     phzy = plot(axh_pmap_y,[-1 sz(2)+1],roi_maxz*[1 1],'-');
%     phxy = plot(axh_pmap_y,roi_maxx*[1 1],[-1 sz(3)+1],'-');
%     hold(axh_pmap_x,'on')
%     phzx = plot(axh_pmap_x,[-1 sz(1)+1],roi_maxz*[1 1],'-');
%     phyx = plot(axh_pmap_x,roi_maxy*[1 1],[-1 sz(3)+1],'-');
%     set([phzy; phzx; phxy; phxz; phyz; phyx],'Color',[1 1 .999])
% 
%     th_z = text(axh_pmap_z,1,1,num2str(roi_maxz));
%     th_y = text(axh_pmap_y,1,1,num2str(roi_maxy));
%     th_x = text(axh_pmap_x,1,1,num2str(roi_maxx));
%     set([th_z; th_y; th_x],'Color',[1 1 .999],'VerticalAlignment','bottom')
%     th_case_z = text(axh_pmap_z,1,sz(2),case_name);
%     set(th_case_z,'Color',[1 1 .999],'VerticalAlignment','top')
% 
%     colormap(axh_pmap_x,'gray')
%     colormap(axh_pmap_y,'gray')
%     colormap(axh_pmap_z,'gray')
%     % set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',max(I(logical(roi_I)))*clim.s2000)
% %     set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',max(I(roi_I))*clim.s2000)
% 
%     I_maxz = I(:,:,roi_maxz);
% 
%    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',quantile(I_maxz(isfinite(I_maxz)),.99,'all')*clim.s2000)
% %    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',max(I(:,:,roi_maxz),[],'all')*clim.s2000)
% %     set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',max(I(:))*clim.s2000)
% %     set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',2*median(I(I>0))*clim.s2000)
%     axis([axh_pmap_x; axh_pmap_y; axh_pmap_z],'off')
% 
    dist_s = diso_sddelta_dist_s;
    contourpars.Nlevels = 10;
    C = contourc(double(dist_s.x),double(dist_s.y),double(dist_s.w'),contourpars.Nlevels);

%     axh_cont = axes('position',[fov.x/(fov.x+fov.y)+.1 fov.z/(fov.y+fov.z)+.1 fov.y/(fov.x+fov.y)-.15 fov.y/(fov.y+fov.z)-.15]);
    axh_cont = axes('position',[.18 .18 .8 .8]);
    hold(axh_cont,'on')

    hcontour = [];
    count = 1;
    while count < length(C)
        numxy = C(2,count);
        xtemp = C(1,count+(1:numxy));
        ytemp = C(2,count+(1:numxy));
        h = plot(axh_cont,xtemp,ytemp,'k-','LineWidth',.5*opt.mplot.lw);
        hcontour = [hcontour; h];
        count = count + numxy + 1;
    end

    xlabel('size, D_{iso} / m^2s^{-1}')
    xlim = clim.mdiso + .05*abs(diff(clim.mdiso))*[-1 1];
    ylabel('shape, D_\Delta^2')
    ylim = clim.msddelta + .05*abs(diff(clim.msddelta))*[-1 1];
    set(axh_cont,'XLim',xlim,'YLim',ylim)
    axis square
    set(axh_cont,'Box','off','TickDir','out','TickLength',.02*[1 1],'FontSize',opt.mplot.fs,'LineWidth',opt.mplot.lw)

    roi_I =logical(mfs.mask);
    ind = roi_I;

    [I_s0,~] = mdm_nii_read(fullfile(nii_path,'dtd_s0.nii.gz'));
    [I_mdiso,~] = mdm_nii_read(fullfile(nii_path,'dtd_mdiso.nii.gz'));
    [I_msddelta,~] = mdm_nii_read(fullfile(nii_path,'dtd_msddelta.nii.gz'));

    dist_s = dist_2d_discrete2smooth([I_mdiso(ind) I_msddelta(ind) I_s0(ind)],dist_s_in);

    C = contourc(double(dist_s.x),double(dist_s.y),double(dist_s.w'),contourpars.Nlevels);

    hcontour = [];
    count = 1;
    while count < length(C)
        numxy = C(2,count);
        xtemp = C(1,count+(1:numxy));
        ytemp = C(2,count+(1:numxy));
        h = plot(axh_cont,xtemp,ytemp,'b-','LineWidth',.5*opt.mplot.lw);
        hcontour = [hcontour; h];
        count = count + numxy + 1;
    end

%     [I_s0,~] = mdm_nii_read(fullfile(nii_path,'dtd_gamma_s0.nii.gz'));
%     [I_mdiso,~] = mdm_nii_read(fullfile(nii_path,'dtd_gamma_MD.nii.gz')); I_mdiso = I_mdiso*1e-9;
%     [I_msddelta,~] = mdm_nii_read(fullfile(nii_path,'dtd_gamma_nmsdaniso.nii.gz'));
% 
%     dist_s = dist_2d_discrete2smooth([I_mdiso(ind) I_msddelta(ind) I_s0(ind)],dist_s_in);
% 
%     C = contourc(double(dist_s.x),double(dist_s.y),double(dist_s.w'),contourpars.Nlevels);
% 
%     hcontour = [];
%     count = 1;
%     while count < length(C)
%         numxy = C(2,count);
%         xtemp = C(1,count+(1:numxy));
%         ytemp = C(2,count+(1:numxy));
%         h = plot(axh_cont,xtemp,ytemp,'g-','LineWidth',.5*opt.mplot.lw);
%         hcontour = [hcontour; h];
%         count = count + numxy + 1;
%     end
% 
%     [I_s0,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_s0.nii.gz'));
%     [I_mdiso,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_MD.nii.gz')); I_mdiso = I_mdiso*1e-9;
%     [I_msddelta,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_nmsdaniso.nii.gz'));
% 
%     dist_s = dist_2d_discrete2smooth([I_mdiso(ind) I_msddelta(ind) I_s0(ind)],dist_s_in);
% 
%     C = contourc(double(dist_s.x),double(dist_s.y),double(dist_s.w'),contourpars.Nlevels);
% 
%     hcontour = [];
%     count = 1;
%     while count < length(C)
%         numxy = C(2,count);
%         xtemp = C(1,count+(1:numxy));
%         ytemp = C(2,count+(1:numxy));
%         h = plot(axh_cont,xtemp,ytemp,'r-','LineWidth',.5*opt.mplot.lw);
%         hcontour = [hcontour; h];
%         count = count + numxy + 1;
%     end
    
%     legend(axh_cont,{'dtd','dtd mean','cov','gamma'},'FontSize',.8*fs)
%     legend(axh_cont,'boxoff','Location','northeast')

    papersize = 8.3*[1 1];
    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    print(fullfile(fileparts(nii_path),'global_dtd.pdf'),'-loose','-dpdf')
%%

    figure(3), clf
    Npanels = 6;
    papersize = 17*[1 1/Npanels/1.618];
%     papersize = 8.3*[1 1.618];

    left = .02;
    dleft = (1-left)/Npanels;
	width = dleft - .02;
    bottom = .18;
%     dbottom = (1-bottom)/Npanels;
    height = 1 - bottom - .01;
    
    lw = 1;
    fs = 6;
    
    Nbins = 50;

    npanel = 1;
    position = [-.005 -.02 1/Npanels 1.02];
%     position = [left bottom+(Npanels-npanel)*dbottom width height];
    edges = linspace(clim.mdiso(1)/1e-9,clim.mdiso(2)/1e-9,Nbins);
    axh_legend = axes('position',position);
    hold(axh_legend,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_MD.nii.gz'));
    [counts_dtd,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_MD.nii.gz'));
    [counts_covariance,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_gamma_MD.nii.gz'));
    [counts_gamma,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dti_lls_MD.nii.gz'));
    [counts_lls,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dti_euler_MD.nii.gz'));
    [counts_euler,~] = histcounts(histdat(ind),edges);
    hh_legend_dtd = histogram(axh_legend,'BinEdges', edges, 'BinCounts', counts_dtd);
    hh_legend_covariance = histogram(axh_legend,'BinEdges', edges, 'BinCounts', counts_covariance);
    hh_legend_gamma = histogram(axh_legend,'BinEdges', edges, 'BinCounts', counts_gamma);
    hh_legend_lls = histogram(axh_legend,'BinEdges', edges, 'BinCounts', counts_lls);
    hh_legend_euler = histogram(axh_legend,'BinEdges', edges, 'BinCounts', counts_euler);
    legend(axh_legend,{'dtd','cov','gamma','lls','euler'},'FontSize',.8*fs)
    legend(axh_legend,'boxoff')
    legend(axh_legend,'Location','southeast')
    ylim = max([counts_dtd(2:(end-1)) counts_covariance(2:(end-1)) counts_gamma(2:(end-1)) counts_lls(2:(end-1)) counts_euler(2:(end-1))])*[-.1 1.1];
    set(axh_legend,'YLim',[0 1e3*ylim(2)],'XTick',[])
    th_case = text(axh_legend,edges(1),1e3*ylim(2),case_name);
    set(th_case,'FontSize',fs,'VerticalAlignment','top')

    
    npanel = 2;
    position = [left+(npanel-1)*dleft bottom width height];
%     position = [left bottom+(Npanels-npanel)*dbottom width height];
    edges = linspace(clim.mdiso(1)/1e-9,clim.mdiso(2)/1e-9,Nbins);
    axh_MD = axes('position',position);
    hold(axh_MD,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_MD.nii.gz'));
    [counts_dtd,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_MD.nii.gz'));
    [counts_covariance,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_gamma_MD.nii.gz'));
    [counts_gamma,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dti_lls_MD.nii.gz'));
    [counts_lls,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dti_euler_MD.nii.gz'));
    [counts_euler,~] = histcounts(histdat(ind),edges);
    hh_MD_dtd = histogram(axh_MD,'BinEdges', edges, 'BinCounts', counts_dtd);
    hh_MD_covariance = histogram(axh_MD,'BinEdges', edges, 'BinCounts', counts_covariance);
    hh_MD_gamma = histogram(axh_MD,'BinEdges', edges, 'BinCounts', counts_gamma);
    hh_MD_lls = histogram(axh_MD,'BinEdges', edges, 'BinCounts', counts_lls);
    hh_MD_euler = histogram(axh_MD,'BinEdges', edges, 'BinCounts', counts_euler);
    ylim = max([counts_dtd(2:(end-1)) counts_covariance(2:(end-1)) counts_gamma(2:(end-1)) counts_lls(2:(end-1)) counts_euler(2:(end-1))])*[-.1 1.1];
    set(axh_MD,'YLim',ylim)

    npanel = 3;
    position = [left+(npanel-1)*dleft bottom width height];
%     position = [left bottom+(Npanels-npanel)*dbottom width height];
    edges = linspace(0,1,Nbins);
    axh_FA = axes('position',position);
    hold(axh_FA,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_FA.nii.gz'));
    [counts_dtd,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_FA.nii.gz'));
    [counts_covariance,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dti_lls_FA.nii.gz'));
    [counts_lls,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dti_euler_FA.nii.gz'));
    [counts_euler,~] = histcounts(histdat(ind),edges);
    hh_FA_dtd = histogram(axh_FA,'BinEdges', edges, 'BinCounts', counts_dtd);
    hh_FA_covariance = histogram(axh_FA,'BinEdges', edges, 'BinCounts', counts_covariance);
    hh_FA_lls = histogram(axh_FA,'BinEdges', edges, 'BinCounts', counts_lls);
    hh_FA_euler = histogram(axh_FA,'BinEdges', edges, 'BinCounts', counts_euler);
    ylim = max([counts_dtd(2:(end-1)) counts_covariance(2:(end-1)) counts_lls(2:(end-1)) counts_euler(2:(end-1))])*[-.1 1.1];
    set(axh_FA,'YLim',ylim)

    npanel = 4;
    position = [left+(npanel-1)*dleft bottom width height];
%     position = [left bottom+(Npanels-npanel)*dbottom width height];
    edges = linspace(0,1,Nbins);
    axh_nmsdaniso = axes('position',position);
    hold(axh_nmsdaniso,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_nmsdaniso.nii.gz'));
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_msddelta.nii.gz'));
    [counts_dtd,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_nmsdaniso.nii.gz'));
    [counts_covariance,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_gamma_nmsdaniso.nii.gz'));
    [counts_gamma,~] = histcounts(histdat(ind),edges);
    hh_shape_dtd = histogram(axh_nmsdaniso,'BinEdges', edges, 'BinCounts', counts_dtd);
    hh_shape_covariance = histogram(axh_nmsdaniso,'BinEdges', edges, 'BinCounts', counts_gamma);
    hh_shape_gamma = histogram(axh_nmsdaniso,'BinEdges', edges, 'BinCounts', counts_covariance);
    ylim = max([counts_dtd(2:(end-1)) counts_covariance(2:(end-1)) counts_gamma(2:(end-1))])*[-.1 1.1];
    set(axh_nmsdaniso,'YLim',ylim)
    

    npanel = 5;
    position = [left+(npanel-1)*dleft bottom width height];
%     position = [left bottom+(Npanels-npanel)*dbottom width height];
    edges = linspace(0,1,Nbins);
    axh_nvdiso = axes('position',position);
    hold(axh_nvdiso,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_nvdiso.nii.gz'));
    [counts_dtd,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_nvdiso.nii.gz'));
    [counts_covariance,~] = histcounts(histdat(ind),edges);
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_gamma_nvdiso.nii.gz'));
    [counts_gamma,~] = histcounts(histdat(ind),edges);
    hh_vsize_dtd = histogram(axh_nvdiso,'BinEdges', edges, 'BinCounts', counts_dtd);
    hh_vsize_gamma = histogram(axh_nvdiso,'BinEdges', edges, 'BinCounts', counts_gamma);
    hh_vsize_covariance = histogram(axh_nvdiso,'BinEdges', edges, 'BinCounts', counts_covariance);
    ylim = max([counts_dtd(2:(end-1)) counts_covariance(2:(end-1)) counts_gamma(2:(end-1))])*[-.1 1.1];
    set(axh_nvdiso,'YLim',ylim)

    npanel = 6;
    position = [left+(npanel-1)*dleft bottom width height];
%     position = [left bottom+(Npanels-npanel)*dbottom width height];
    edges = linspace(0,1,Nbins);
    axh_fractions = axes('position',position);
    hold(axh_fractions,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_fractions.nii.gz'));
    histdat1 = squeeze(double(histdat(1,:,:,:))./sum(histdat));
    histdat2 = squeeze(double(histdat(2,:,:,:))./sum(histdat));
    histdat3 = squeeze(double(histdat(3,:,:,:))./sum(histdat));
    [counts1,~] = histcounts(histdat1(ind),edges);
    [counts2,~] = histcounts(histdat2(ind),edges);
    [counts3,~] = histcounts(histdat3(ind),edges);
    hh_fraction3 = histogram(axh_fractions,'BinEdges', edges, 'BinCounts', counts3);
    hh_fraction2 = histogram(axh_fractions,'BinEdges', edges, 'BinCounts', counts2);
    hh_fraction1 = histogram(axh_fractions,'BinEdges', edges, 'BinCounts', counts1);
    ylim = max([counts1(2:(end-1)) counts2(2:(end-1)) counts3(2:(end-1))])*[-.1 1.1];
    set(axh_fractions,'YLim',ylim)

    
    set([hh_legend_dtd; hh_legend_covariance; hh_legend_gamma; hh_legend_lls; hh_legend_euler;...
        hh_MD_dtd; hh_MD_covariance; hh_MD_gamma; hh_MD_lls; hh_MD_euler;...
        hh_shape_dtd; hh_shape_covariance; hh_shape_gamma;...
        hh_vsize_dtd; hh_vsize_covariance; hh_vsize_gamma;...
        hh_FA_dtd; hh_FA_covariance; hh_FA_lls; hh_FA_euler;...
        hh_fraction1; hh_fraction2; hh_fraction3],'DisplayStyle','stairs')
    set([hh_legend_dtd; hh_MD_dtd; hh_shape_dtd; hh_vsize_dtd; hh_FA_dtd; hh_fraction3],'EdgeColor','b','LineWidth',1*lw)
    set([hh_legend_gamma; hh_MD_gamma; hh_shape_gamma; hh_vsize_gamma; hh_fraction2],'EdgeColor','g','LineWidth',.75*lw)
    set([hh_legend_covariance; hh_MD_covariance; hh_shape_covariance; hh_vsize_covariance; hh_FA_covariance; hh_fraction1],'EdgeColor','r','LineWidth',.5*lw)
    set([hh_legend_lls; hh_MD_lls; hh_FA_lls],'EdgeColor','c')
    set([hh_legend_euler; hh_MD_euler; hh_FA_euler],'EdgeColor','m')
    
    set([axh_legend; axh_MD; axh_nmsdaniso; axh_nvdiso; axh_FA; axh_fractions],'Box','off','TickDir','out',...
        'LineWidth',lw,'FontSize',fs,'YTick',[])
    xlabel(axh_MD,'mean size, MD / 10^{-9} m^2s^{-1}')
    xlabel(axh_FA,'shape x alignment, FA')
    xlabel(axh_nmsdaniso,'mean-square shape  <D_{aniso}^2>/<D_{iso}>^2')
    xlabel(axh_nvdiso,'variance of sizes, Var(D_{iso})/<D_{iso}>^2')
    xlabel(axh_fractions,'dtd fractions thin (R), thick (G), big (B)')
    
    
    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    print(fullfile(fileparts(nii_path),['global_histograms.pdf']),'-loose','-dpdf')
    
    %%
%     figure(10), clf
%     wm_mask = histdat1 > .5;
%     nk = 13;
% %     imagesc(wm_mask(:,:,nk)')
%     imagesc(histdat1(:,:,nk)')

end

