clear all

data_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419';
data_dir = dir(fullfile(data_path,'*'));
rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois/DT';
% rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois/Yuan20200320';
% rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois/Yuan20200306';
% rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois_20200304';
% % rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois_20200215';
% rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois/Yuan20200204';

% data_path = '/Users/daniel/Dropbox/NMRdata/GE_MSKCC_breast/processing_dicom_manmask/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'20*'));
% rois_path = '/Users/daniel/Dropbox/NMRdata/GE_MSKCC_breast/processing_dicom_manmask/rois_DT';

% data_path = '/Users/daniel/Dropbox/NMRdata/GE_Spectrum_1/processing_dicom_manmask/acquisition_date_time';
% data_dir = dir(fullfile(data_path,'20*'));
% rois_path = '/Users/daniel/Dropbox/NMRdata/GE_Spectrum_1/processing_dicom_manmask/rois_DT';

% Make paths to nii_xps, boostraps, and maps folders
method = 'dtd';

[~,rois_name,~] = fileparts(rois_path);
data_names = cell(size(data_dir,1),1);
for ndata = 1:size(data_dir,1)
    data_names{ndata} = data_dir(ndata).name;
end
data_names = data_names(~contains(data_names,'.'));
% data_names = data_names(71);
% data_names = data_names(contains(data_names,{'P0338456'}));

nii_fns = cell(0,0);
pmaps_paths = cell(0,0);
bs_paths = cell(0,0);
for ndata = 1:numel(data_names)
    nii_fns{1+numel(nii_fns)} = fullfile(data_path,data_names{ndata},'rwi','nii_xps','data_mc.nii.gz');
    rwi_path = fileparts(fileparts(nii_fns{ndata}));
    pmaps_paths{1+numel(pmaps_paths)} = fullfile(rwi_path,'pmaps');
    bs_paths{1+numel(bs_paths)} = fullfile(fileparts(fileparts(nii_fns{ndata})),method,'bootstraps');        
end

[roi_plots_path,roi_plots_name,~] = fileparts(rois_path);
roi_plots_path = fullfile(fileparts(roi_plots_path),'roi_plots',roi_plots_name);

pmap_ext = '.nii.gz';
roi_prefix = 'lesion';
% roi_prefix = 'wm';

lw = 1;
fs = 8;

% Define color limits for parameter maps
clim.s0 = 1*[0 1]; % Multiplied with s0 below
clim.s2000 = 1*[0 1]; % Multiplied with s0 below
clim.mdiso = 3.5e-9*[0 1]; %clim.mdiso = 1e-9*[0 1];
clim.msddelta = 1*[0 1];
clim.vdiso = .3*3e-9^2*[0 1]; %clim.vdiso = .3*1e-9^2*[0 1];
clim.vsddelta = .15*[0 1];
clim.cvdisosddelta = .1*3e-9*1*[-1 1]; %clim.cvdisosddelta = .2*1e-9*1*[-1 1];
Nbins = 50;

for ndata = 1:numel(nii_fns)
    
    data_name = data_names{ndata};
    nii_fn = nii_fns{ndata};
    roi_path = fullfile(rois_path,data_name);
    
    bs_path = bs_paths{ndata};
    pmaps_path = pmaps_paths{ndata};    
        
    
    % Delete old version of ROI pdf files
%     roi_pdf_dir = dir(fullfile(fileparts(pmaps_path),[roi_prefix '*.pdf']));
%     for npdf = 1:numel(roi_pdf_dir)
%         delete(fullfile(roi_pdf_dir(npdf).folder,roi_pdf_dir(npdf).name))
%     end
    
%     roi_dir = dir(fullfile(roi_path,[roi_prefix '*Roi*']));
    roi_dir = dir(fullfile(roi_path,[roi_prefix '*']));    
    
    if exist(pmaps_path,'dir') ~= 7, disp([data_name ' pmaps_path missing']), continue, end
    if exist(bs_path,'dir') ~= 7, disp([data_name ' bs_path missing']), continue, end
    if exist(nii_fn,'file') ~= 2, disp([data_name ' nii_fn missing']), continue, end
    if isempty(roi_dir), disp([data_name ' roi_dir missing']), continue, end

    % Extract and clean up roi_name for later use in output file names
    roi_nameext = roi_dir.name;
    roi_fn = fullfile(roi_path,roi_nameext);
    [roi_path,roi_name,roi_ext] = fileparts(roi_fn);
    ind_dot = find(roi_name == '.');
    if ~isempty(ind_dot)
        roi_name = roi_name(1:(ind_dot-1));
    end
    pattern = '_Roi'; % From Wuhan naming convention
    ind_pattern = strfind(roi_name,pattern);
    if ~isempty(ind_pattern)
        roi_name = roi_name(1:(ind_pattern-1));
    end
    
    if exist(roi_fn,'file') ~= 2, disp([data_name ' roi_fn missing']), continue, end
    
    msf_mkdir(roi_plots_path)
    fig_fn = fullfile(roi_plots_path,[data_name '_' roi_name '_b0_b2000_MD_FA.pdf']);
    
    opt = mdm_opt();
    opt = mplot_opt(opt);
    opt.roiplot.data_name = data_name;

    [roi_I,roi_h] = mdm_nii_read(roi_fn);
    roi_I = logical(roi_I);

    roi_projx = squeeze(sum(sum(roi_I,2),3));
    roi_maxx = max(find(roi_projx==max(roi_projx)));
    roi2d_x = squeeze(roi_I(roi_maxx,:,:));

    roi_projy = squeeze(sum(sum(roi_I,1),3));
    roi_maxy = max(find(roi_projy==max(roi_projy)));
    roi2d_y = squeeze(roi_I(:,roi_maxy,:));

    roi_projz = squeeze(sum(sum(roi_I,1),2));
    roi_maxz = max(find(roi_projz==max(roi_projz)));
    roi2d_z = squeeze(roi_I(:,:,roi_maxz));
  

    pmap_name = 'dtd_gamma_s2000';
    pmap_fn = fullfile(pmaps_path,[pmap_name pmap_ext]);
    
    [I,nii_h] = mdm_nii_read(pmap_fn);
    I = double(I);
    sz = size(I);
    
    figure(2), clf
    
    fov.x = nii_h.dim(2)*nii_h.pixdim(2);
    fov.y = nii_h.dim(3)*nii_h.pixdim(3);
    fov.z = nii_h.dim(4)*nii_h.pixdim(4);

    position_hist_0   = [fov.x/(fov.x+fov.y)+.05 .14 fov.y/(fov.x+fov.y)-.1 fov.y/(fov.y+fov.z)-.14];
    position_pmap_x_0 = [fov.x/(fov.x+fov.y) fov.y/(fov.y+fov.z) fov.y/(fov.x+fov.y) fov.z/(fov.y+fov.z)];
    position_pmap_y_0 = [0 fov.y/(fov.y+fov.z) fov.x/(fov.x+fov.y) fov.z/(fov.y+fov.z)];
    position_pmap_z_0 = [0 0 fov.x/(fov.x+fov.y) fov.y/(fov.y+fov.z)];

    Npmaps = 4;
    npmap = 2;
    
    position_hist     = position_hist_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_x   = position_pmap_x_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_y   = position_pmap_y_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_z   = position_pmap_z_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];

    
    im2d_x = squeeze(I(roi_maxx,:,:));
    im2d_y = squeeze(I(:,roi_maxy,:));
    im2d_z = squeeze(I(:,:,roi_maxz));

    axh_pmap_x = axes('position',position_pmap_x);
    axh_pmap_y = axes('position',position_pmap_y);
    axh_pmap_z = axes('position',position_pmap_z);

    imagesc(axh_pmap_x,im2d_x')
    imagesc(axh_pmap_y,im2d_y')
    imagesc(axh_pmap_z,im2d_z')
    
    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'YDir','normal')

    colormap(axh_pmap_x,'gray')
    colormap(axh_pmap_y,'gray')
    colormap(axh_pmap_z,'gray')

    plot_roi(roi2d_x',[1 0 0], .5*lw, 1, axh_pmap_x);
    plot_roi(roi2d_y',[1 0 0], .5*lw, 1, axh_pmap_y);
    plot_roi(roi2d_z',[1 0 0], .5*lw, 1, axh_pmap_z);

    I_maxz = I(:,:,roi_maxz);

    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',quantile(I_maxz(I_maxz>0),.999,'all')*clim.s2000)
    axis([axh_pmap_x; axh_pmap_y; axh_pmap_z],'off')

        
    axh_hist = axes('position',position_hist);
    hold(axh_hist,'on')
    [histdat,~] = mdm_nii_read(pmap_fn);
    edges = linspace(0,max(histdat(roi_I),[],'all'),Nbins);
    [counts,~] = histcounts(histdat(roi_I),edges);
    hh = histogram(axh_hist,'BinEdges', edges, 'BinCounts', counts);
    ylim = max(counts(2:(end-1)))*[-.1 1.1];
    set(axh_hist,'YLim',ylim)
    
    set(hh,'DisplayStyle','stairs','EdgeColor','b')
    set(axh_hist,'Box','off','TickDir','out','TickLength',.02*[1 1],'LineWidth',lw,'FontSize',fs,'YTick',[])

    set([axh_hist],'Box','off','TickDir','out','TickLength',.02*[1 1],'LineWidth',lw,'FontSize',fs,'YTick',[])
    xlabel(axh_hist,'S(b2000)','FontSize',fs)
    

    pmap_name = 'dti_euler_s0';    
    pmap_fn = fullfile(pmaps_path,[pmap_name pmap_ext]);
    
    npmap = 1;
    
    position_hist     = position_hist_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_x   = position_pmap_x_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_y   = position_pmap_y_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_z   = position_pmap_z_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];

    [I,nii_h] = mdm_nii_read(pmap_fn);
    I = double(I);
    sz = size(I);
    
    
    fov.x = nii_h.dim(2)*nii_h.pixdim(2);
    fov.y = nii_h.dim(3)*nii_h.pixdim(3);
    fov.z = nii_h.dim(4)*nii_h.pixdim(4);

    im2d_x = squeeze(I(roi_maxx,:,:));
    im2d_y = squeeze(I(:,roi_maxy,:));
    im2d_z = squeeze(I(:,:,roi_maxz));

    axh_pmap_x = axes('position',position_pmap_x);
    axh_pmap_y = axes('position',position_pmap_y);
    axh_pmap_z = axes('position',position_pmap_z);

    imagesc(axh_pmap_x,im2d_x')
    imagesc(axh_pmap_y,im2d_y')
    imagesc(axh_pmap_z,im2d_z')
    
    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'YDir','normal')

    hold(axh_pmap_z,'on')
    phxz = plot(axh_pmap_z,roi_maxx*[1 1],[-1 sz(1)+1],'-');
    phyz = plot(axh_pmap_z,[-1 sz(2)+1],roi_maxy*[1 1],'-');
    hold(axh_pmap_y,'on')
    phzy = plot(axh_pmap_y,[-1 sz(2)+1],roi_maxz*[1 1],'-');
    phxy = plot(axh_pmap_y,roi_maxx*[1 1],[-1 sz(3)+1],'-');
    hold(axh_pmap_x,'on')
    phzx = plot(axh_pmap_x,[-1 sz(1)+1],roi_maxz*[1 1],'-');
    phyx = plot(axh_pmap_x,roi_maxy*[1 1],[-1 sz(3)+1],'-');
    set([phzy; phzx; phxy; phxz; phyz; phyx],'Color',[1 1 .999])
    
    th_x = text(axh_pmap_x,.5,sz(3)+.5,['x = ' num2str(roi_maxx)]);
    th_y = text(axh_pmap_y,.5,sz(3)+.5,['y = ' num2str(roi_maxy)]);
    th_z = text(axh_pmap_z,.5,sz(2)+.5,['z = ' num2str(roi_maxz)]);
    set([th_z; th_y; th_x],'Color',[1 1 .999],'VerticalAlignment','top')
    th_data_z = text(axh_pmap_z,.5,.5,data_name);
    set(th_data_z,'Color',[1 1 .999],'VerticalAlignment','bottom')
    set([th_z; th_y; th_x; th_data_z],'Color',[1 1 .999],'FontSize',fs)

    colormap(axh_pmap_x,'gray')
    colormap(axh_pmap_y,'gray')
    colormap(axh_pmap_z,'gray')

    plot_roi(roi2d_x',[1 0 0], .5*lw, 1, axh_pmap_x);
    plot_roi(roi2d_y',[1 0 0], .5*lw, 1, axh_pmap_y);
    plot_roi(roi2d_z',[1 0 0], .5*lw, 1, axh_pmap_z);

    I_maxz = I(:,:,roi_maxz);

    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',quantile(I_maxz(I_maxz>0),.999,'all')*clim.s0)
%     set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',clim_temp)
    axis([axh_pmap_x; axh_pmap_y; axh_pmap_z],'off')    
            
    axh_hist = axes('position',position_hist);
    hold(axh_hist,'on')
    [histdat,~] = mdm_nii_read(pmap_fn);
    edges = linspace(0,max(histdat(roi_I),[],'all'),Nbins);
    [counts,~] = histcounts(histdat(roi_I),edges);
    hh = histogram(axh_hist,'BinEdges', edges, 'BinCounts', counts);
    ylim = max(counts(2:(end-1)))*[-.1 1.1];
    set(axh_hist,'YLim',ylim)
    
    set(hh,'DisplayStyle','stairs','EdgeColor','b')
    set(axh_hist,'Box','off','TickDir','out','TickLength',.02*[1 1],'LineWidth',lw,'FontSize',fs,'YTick',[])

    xlabel(axh_hist,'S(b0)','FontSize',fs)
    
    
    %%   
    pmap_name = 'dti_euler_MD';
    clim_temp = clim.mdiso*1e9;
    pmap_fn = fullfile(pmaps_path,[pmap_name pmap_ext]);
    
    npmap = 3;
    
    position_hist     = position_hist_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_x   = position_pmap_x_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_y   = position_pmap_y_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_z   = position_pmap_z_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    
    [I,nii_h] = mdm_nii_read(pmap_fn);
    I = double(I);
    sz = size(I);
    
    
    fov.x = nii_h.dim(2)*nii_h.pixdim(2);
    fov.y = nii_h.dim(3)*nii_h.pixdim(3);
    fov.z = nii_h.dim(4)*nii_h.pixdim(4);

    im2d_x = squeeze(I(roi_maxx,:,:));
    im2d_y = squeeze(I(:,roi_maxy,:));
    im2d_z = squeeze(I(:,:,roi_maxz));

    axh_pmap_x = axes('position',position_pmap_x);
    axh_pmap_y = axes('position',position_pmap_y);
    axh_pmap_z = axes('position',position_pmap_z);

    imagesc(axh_pmap_x,im2d_x')
    imagesc(axh_pmap_y,im2d_y')
    imagesc(axh_pmap_z,im2d_z')
    
    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'YDir','normal')

    colormap(axh_pmap_x,'gray')
    colormap(axh_pmap_y,'gray')
    colormap(axh_pmap_z,'gray')

    plot_roi(roi2d_x',[1 0 0], .5*lw, 1, axh_pmap_x);
    plot_roi(roi2d_y',[1 0 0], .5*lw, 1, axh_pmap_y);
    plot_roi(roi2d_z',[1 0 0], .5*lw, 1, axh_pmap_z);

    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',clim_temp)
    axis([axh_pmap_x; axh_pmap_y; axh_pmap_z],'off')    
            
    axh_hist = axes('position',position_hist);
    hold(axh_hist,'on')
    [histdat,~] = mdm_nii_read(pmap_fn);
    edges = linspace(clim_temp(1),clim_temp(2),Nbins);
    [counts,~] = histcounts(histdat(roi_I),edges);
    hh = histogram(axh_hist,'BinEdges', edges, 'BinCounts', counts);
    ylim = max(counts(2:(end-1)))*[-.1 1.1];
    set(axh_hist,'XLim',clim_temp,'YLim',ylim)
    
    set(hh,'DisplayStyle','stairs','EdgeColor','b')
    set(axh_hist,'Box','off','TickDir','out','TickLength',.02*[1 1],'LineWidth',lw,'FontSize',fs,'YTick',[])

    xlabel(axh_hist,'MD / 10^{-9}m^2s^{-1}','FontSize',fs)
    
    
    %%   
    pmap_name = 'dti_euler_FA_u_rgb';
    clim_temp = [0 1];
    pmap_fn = fullfile(pmaps_path,[pmap_name pmap_ext]);
    
    npmap = 4;
    
    position_hist     = position_hist_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_x   = position_pmap_x_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_y   = position_pmap_y_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    position_pmap_z   = position_pmap_z_0.*[1/Npmaps 1 1/Npmaps 1] + [(npmap-1)/Npmaps 0 0 0];
    
    [I,nii_h] = mdm_nii_read(pmap_fn);
    I = double(I);
    sz = size(I);
    
    sz = sz(2:4);
    I_r = reshape(I(1,:,:,:),sz)/255;
    I_g = reshape(I(2,:,:,:),sz)/255;
    I_b = reshape(I(3,:,:,:),sz)/255;
    
    
    fov.x = nii_h.dim(2)*nii_h.pixdim(2);
    fov.y = nii_h.dim(3)*nii_h.pixdim(3);
    fov.z = nii_h.dim(4)*nii_h.pixdim(4);

    im2d_x = cat(3,squeeze(I_r(roi_maxx,:,:))',squeeze(I_g(roi_maxx,:,:))',squeeze(I_b(roi_maxx,:,:))');
    im2d_y = cat(3,squeeze(I_r(:,roi_maxy,:))',squeeze(I_g(:,roi_maxy,:))',squeeze(I_b(:,roi_maxy,:))');
    im2d_z = cat(3,squeeze(I_r(:,:,roi_maxz))',squeeze(I_g(:,:,roi_maxz))',squeeze(I_b(:,:,roi_maxz))');
    
    axh_pmap_x = axes('position',position_pmap_x);
    axh_pmap_y = axes('position',position_pmap_y);
    axh_pmap_z = axes('position',position_pmap_z);

    imagesc(axh_pmap_x,im2d_x)
    imagesc(axh_pmap_y,im2d_y)
    imagesc(axh_pmap_z,im2d_z)
    
    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'YDir','normal')

    plot_roi(roi2d_x',[1 1 .999], .5*lw, 1, axh_pmap_x);
    plot_roi(roi2d_y',[1 1 .999], .5*lw, 1, axh_pmap_y);
    plot_roi(roi2d_z',[1 1 .999], .5*lw, 1, axh_pmap_z);

    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',clim_temp)
    axis([axh_pmap_x; axh_pmap_y; axh_pmap_z],'off')
    
    pmap_name = 'dti_euler_FA';
    pmap_fn = fullfile(pmaps_path,[pmap_name pmap_ext]);
    
    position = [fov.x/(fov.x+fov.y)+.05 .12 fov.y/(fov.x+fov.y)-.1 fov.y/(fov.y+fov.z)-.13];
            
    axh_hist = axes('position',position_hist);
    hold(axh_hist,'on')
    [histdat,~] = mdm_nii_read(pmap_fn);
    if strcmp(class(histdat),'uint8')
        histdat = double(histdat)/255;
    end
    edges = linspace(clim_temp(1),clim_temp(2),Nbins);
    [counts,~] = histcounts(histdat(roi_I),edges);
    hh = histogram(axh_hist,'BinEdges', edges, 'BinCounts', counts);
    ylim = max(counts(2:(end-1)))*[-.1 1.1];
    set(axh_hist,'XLim',clim_temp,'YLim',ylim)
    
    set(hh,'DisplayStyle','stairs','EdgeColor','b')
    set(axh_hist,'Box','off','TickDir','out','TickLength',.02*[1 1],'LineWidth',lw,'FontSize',fs,'YTick',[])

    xlabel(axh_hist,'FA','FontSize',fs)
    
    papersize = 8.3*[1*Npmaps (fov.y+fov.z)/(fov.y+fov.x)];
    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    print(fig_fn,'-loose','-dpdf') 

end

