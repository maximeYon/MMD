clear all


nii_paths = cell(0);
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/data_processing_by_rwi/RWI/RWI 10 - RWI10/Mr Breast Research Exam/MDD intermediate FS - 9/res_dtd/MDD_intermediate_FS/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/data_processing_by_rwi/RWI/RWI 9 - RWI9/Mr Breast Research Exam/MDD intermediate FS - 4/res_dtd/MDD_intermediate_FS/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/data_processing_by_rwi/RWI/RWI 8 - RWI8/Mr Breast Research Exam/MDD intermediate FS - 9/res_dtd/MDD_intermediate_FS/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/data_processing_by_rwi/RWI/RWI 7 - RWI7/Mr Breast Research Exam/MDD intermediate FS - 8/res_dtd/MDD_intermediate_FS/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study29718/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30297/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30307/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30461/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_FS_Study30498/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_intermediate_SpFS_Study30297/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_long_Study29718/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_long_Study30498/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/MSKCC/Patients/mdd/MDD_long_Study30512/dtd/maps/nii';
nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/United/data4dtdpaper/mdd/ONION_MEN_297359/dtd/maps/nii';
nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/United/data4dtdpaper/mdd/p0252374_145532/dtd/maps/nii';
nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/United/data4dtdpaper/mdd/p0278473_102434/dtd/maps/nii';
nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/United/data4dtdpaper/mdd/p0288182_143537/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/GE/2018-11-22_Premier_corrected_waveforms/mdd/long/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/GE/2018-11-22_Premier_corrected_waveforms/mdd/intermediate/dtd/maps/nii';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0240306/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0275752/rwi/pmaps';
% nii_paths{1+numel(nii_paths)} = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419/P0333068/rwi/pmaps';

nii_ext = '.nii.gz';
nii_name = 'dtd_s2000';
roi_name = 'roi_lesion';
% roi_name = 'roi_wm';
% roi_name = 'roi_csf';
% roi_name = 'roi_gm';
% roi_name = 'roi_tha';
% roi_name = 'roi_cb';

for nnii = 1:numel(nii_paths)
    nii_path = nii_paths{nnii};
    [~,data_name,~] = fileparts(fileparts(fileparts(fileparts(nii_path))));

    bs_path = fullfile(fileparts(fileparts(nii_path)),'bootstraps');
%     bs_path = fullfile(fileparts(nii_path),'dtd','bootstraps');

    nii_fn = fullfile(nii_path,[nii_name nii_ext]);
    roi_fn = fullfile(nii_path,[roi_name nii_ext]);

    opt = mdm_opt();
    opt = mplot_opt(opt);

    [I,nii_h] = mdm_nii_read(nii_fn);
    [roi_I,roi_h] = mdm_nii_read(roi_fn);
    roi_I = logical(roi_I);


    %Plot distributions for ROI

    clim.mdiso = [0 3.5]*1e-9;
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
            [dpar,dperp,~,~,w] = dtd_4d_m2pars(mfs.m);
            [diso,~,~,~,~,sddelta] = dtd_pars2dpars(dpar,dperp); 

            nn = size(w,4); 
            ind = logical(repmat(roi_I,[1 1 1 nn]));
            dist_s = dist_2d_discrete2smooth([diso(ind) sddelta(ind) w(ind)],dist_s_in);
            break
        end
    end

    bs_w = NaN*ones([contourpars.Nx,contourpars.Ny,numel(bsno)]);
    tic
    parfor nbs = 1:numel(bsno)
        mfs_fn   = fullfile(bs_path,num2str(bsno(nbs)),'mfs.mat');
        if exist(mfs_fn,'file')==2

            mfs = mdm_mfs_load(mfs_fn);
            [dpar,dperp,~,~,w] = dtd_4d_m2pars(mfs.m);
            [diso,~,~,~,~,sddelta] = dtd_pars2dpars(dpar,dperp); 

            dist_s_temp = dist_2d_discrete2smooth([diso(ind) sddelta(ind) w(ind)],dist_s_in);        
            bs_w(:,:,nbs) = dist_s_temp.w;                
        end
    end
    toc

    dist_s.w = nanmedian(bs_w,3);
    diso_sddelta_dist_s = dist_s;
    clear bs_w
    %%
    
    figure(2), clf

    clim.s2000 = .5*[0 1];
    fov.x = nii_h.dim(2)*nii_h.pixdim(2);
    fov.y = nii_h.dim(3)*nii_h.pixdim(3);
    fov.z = nii_h.dim(4)*nii_h.pixdim(4);

    roi_projx = squeeze(sum(sum(roi_I,2),3));
    roi_maxx = max(find(roi_projx==max(roi_projx)));
    im2d_x = squeeze(I(roi_maxx,:,:));
    roi2d_x = squeeze(roi_I(roi_maxx,:,:));

    roi_projy = squeeze(sum(sum(roi_I,1),3));
    roi_maxy = max(find(roi_projy==max(roi_projy)));
    im2d_y = squeeze(I(:,roi_maxy,:));
    roi2d_y = squeeze(roi_I(:,roi_maxy,:));

    roi_projz = squeeze(sum(sum(roi_I,1),2));
    roi_maxz = max(find(roi_projz==max(roi_projz)));
    im2d_z = squeeze(I(:,:,roi_maxz));
    roi2d_z = squeeze(roi_I(:,:,roi_maxz));

    axh_pmap_z = axes('position',[0 fov.z/(fov.y+fov.z) fov.x/(fov.x+fov.y) fov.y/(fov.y+fov.z)]);
    imagesc(axh_pmap_z,im2d_z')
    set(axh_pmap_z,'YDir','normal')
    plot_roi(roi2d_z');

    axh_pmap_x = axes('position',[fov.x/(fov.x+fov.y) 0 fov.y/(fov.x+fov.y) fov.z/(fov.y+fov.z)]);
    imagesc(axh_pmap_x,im2d_x')
    set(axh_pmap_x,'YDir','normal')
    plot_roi(roi2d_x');

    axh_pmap_y = axes('position',[0 0 fov.x/(fov.x+fov.y) fov.z/(fov.y+fov.z)]);
    imagesc(axh_pmap_y,im2d_y')
    set(axh_pmap_y,'YDir','normal')
    plot_roi(roi2d_y');

    colormap(axh_pmap_x,'gray')
    colormap(axh_pmap_y,'gray')
    colormap(axh_pmap_z,'gray')
    % set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',max(I(logical(roi_I)))*clim.s2000)
    set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',max(I(:))*clim.s2000)
%     set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',2*median(I(I>0))*clim.s2000)
    axis([axh_pmap_x; axh_pmap_y; axh_pmap_z],'off')

    dist_s = diso_sddelta_dist_s;
    contourpars.Nlevels = 10;
    C = contourc(double(dist_s.x),double(dist_s.y),double(dist_s.w'),contourpars.Nlevels);

    axh_cont = axes('position',[fov.x/(fov.x+fov.y)+.1 fov.z/(fov.y+fov.z)+.1 fov.y/(fov.x+fov.y)-.15 fov.y/(fov.y+fov.z)-.15]);
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

    papersize = 17.5*[1 (fov.y+fov.z)/(fov.y+fov.x)];
    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    print(fullfile(fileparts(nii_path),[data_name '_' roi_name '_s2000']),'-loose','-dpdf')
%%

    figure(3), clf
    Npanels = 5;
    left = .1;
	width = .85;
    bottom = .08;
    dbottom = (1-bottom)/Npanels;
    height = dbottom - .08;
    
    lw = 1;
    fs = 8;
    
    Nbins = 50;

    npanel = 1;
    edges = linspace(clim.mdiso(1)/1e-9,clim.mdiso(2)/1e-9,Nbins);
    axh_MD = axes('position',[left bottom+(Npanels-npanel)*dbottom width height]);
    hold(axh_MD,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_MD.nii.gz'));
    hh_MD_dtd = histogram(axh_MD,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_MD.nii.gz'));
%     hh_MD_covariance = histogram(axh_MD,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_gamma_MD.nii.gz'));
%     hh_MD_gamma = histogram(axh_MD,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dti_lls_MD.nii.gz'));
%     hh_MD_lls = histogram(axh_MD,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dti_euler_MD.nii.gz'));
%     hh_MD_euler = histogram(axh_MD,histdat(ind),edges);

    npanel = 2;
    edges = linspace(0,1,Nbins);
    axh_FA = axes('position',[left bottom+(Npanels-npanel)*dbottom width height]);
    hold(axh_FA,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_FA.nii.gz'));
    hh_FA_dtd = histogram(axh_FA,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_FA.nii.gz'));
%     hh_FA_covariance = histogram(axh_FA,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dti_lls_FA.nii.gz'));
%     hh_FA_lls = histogram(axh_FA,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dti_euler_FA.nii.gz'));
%     hh_FA_euler = histogram(axh_FA,histdat(ind),edges);

    npanel = 3;
    edges = linspace(0,1,Nbins);
    axh_nmsdaniso = axes('position',[left bottom+(Npanels-npanel)*dbottom width height]);
    hold(axh_nmsdaniso,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_nmsdaniso.nii.gz'));
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_msddelta.nii.gz'));
    hh_shape_dtd = histogram(axh_nmsdaniso,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_nmsdaniso.nii.gz'));
%     hh_shape_covariance = histogram(axh_nmsdaniso,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_gamma_nmsdaniso.nii.gz'));
%     hh_shape_gamma = histogram(axh_nmsdaniso,histdat(ind),edges);
    

    npanel = 4;
    edges = linspace(0,1,Nbins);
    axh_nvdiso = axes('position',[left bottom+(Npanels-npanel)*dbottom width height]);
    hold(axh_nvdiso,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_nvdiso.nii.gz'));
    hh_vsize_dtd = histogram(axh_nvdiso,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_covariance_nvdiso.nii.gz'));
%     hh_vsize_covariance = histogram(axh_nvdiso,histdat(ind),edges);
%     [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_gamma_nvdiso.nii.gz'));
%     hh_vsize_gamma = histogram(axh_nvdiso,histdat(ind),edges);

    npanel = 5;
    edges = linspace(0,1,Nbins);
    axh_fractions = axes('position',[left bottom+(Npanels-npanel)*dbottom width height]);
    hold(axh_fractions,'on')
    [histdat,~] = mdm_nii_read(fullfile(nii_path,'dtd_fractions.nii.gz'));
    histdat1 = squeeze(double(histdat(1,:,:,:))./sum(histdat));
    histdat2 = squeeze(double(histdat(2,:,:,:))./sum(histdat));
    histdat3 = squeeze(double(histdat(3,:,:,:))./sum(histdat));
    hh_fraction1 = histogram(axh_fractions,histdat1(ind),edges);
    hh_fraction2 = histogram(axh_fractions,histdat2(ind),edges);
    hh_fraction3 = histogram(axh_fractions,histdat3(ind),edges);

    
%     set([hh_MD_dtd; hh_MD_covariance; hh_MD_gamma; hh_MD_lls; hh_MD_euler;...
%         hh_shape_dtd; hh_shape_covariance; hh_shape_gamma;...
%         hh_vsize_dtd; hh_vsize_covariance; hh_vsize_gamma;...
%         hh_FA_dtd; hh_FA_covariance; hh_FA_lls; hh_FA_euler;...
%         hh_fraction1; hh_fraction2; hh_fraction3],'DisplayStyle','stairs')
%     set([hh_MD_dtd; hh_shape_dtd; hh_vsize_dtd; hh_FA_dtd; hh_fraction3],'EdgeColor','b')
%     set([hh_MD_covariance; hh_shape_covariance; hh_vsize_covariance; hh_FA_covariance; hh_fraction1],'EdgeColor','r')
%     set([hh_MD_gamma; hh_shape_gamma; hh_vsize_gamma; hh_fraction2],'EdgeColor','g')
%     set([hh_MD_lls; hh_FA_lls],'EdgeColor','c')
%     set([hh_MD_euler; hh_FA_euler],'EdgeColor','m')
%     
%     set([axh_MD; axh_nmsdaniso; axh_nvdiso; axh_FA; axh_fractions],'Box','off','TickDir','out',...
%         'LineWidth',lw,'FontSize',fs)
    set([hh_MD_dtd; hh_shape_dtd; hh_vsize_dtd; hh_FA_dtd;...
        hh_fraction1; hh_fraction2; hh_fraction3],'DisplayStyle','stairs')
    set([hh_MD_dtd; hh_shape_dtd; hh_vsize_dtd; hh_FA_dtd; hh_fraction3],'EdgeColor','b')
    set([hh_fraction1],'EdgeColor','r')
    set([hh_fraction2],'EdgeColor','g')
    
    set([axh_MD; axh_nmsdaniso; axh_nvdiso; axh_FA; axh_fractions],'Box','off','TickDir','out',...
        'LineWidth',lw,'FontSize',fs)
    xlabel(axh_MD,'mean size, MD / 10^{-9} m^2s^{-1}')
    xlabel(axh_FA,'shape x alignmment, FA')
    xlabel(axh_nmsdaniso,'normalized mean-square shape  <D_{aniso}^2>/<D_{iso}>^2')
    xlabel(axh_nvdiso,'normalized variance of sizes, Var(D_{iso})/<D_{iso}>^2')
    xlabel(axh_fractions,'fractions thin (R), thick (G), big (B)')
    papersize = 8.3*[1 1.618];
    set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
    print(fullfile(fileparts(nii_path),[data_name '_' roi_name '_histograms']),'-loose','-dpdf')

end

return
axh_scatter = axes('position',[0.6 0.3 0.35 0.4]);
plot3(I_mdiso(logical(roi_I))/1e-9,I_nvdiso(logical(roi_I)),I_msddelta(logical(roi_I)),'.')
axis([0 3 0 1 0 1])
axis square
hold on
plot3(I_mdiso(logical(roi_I))/1e-9,0*I_nvdiso(logical(roi_I))+1,I_msddelta(logical(roi_I)),'k.')
plot3(0*I_mdiso(logical(roi_I))/1e-9,I_nvdiso(logical(roi_I)),I_msddelta(logical(roi_I)),'k.')
plot3(I_mdiso(logical(roi_I))/1e-9,I_nvdiso(logical(roi_I)),0*I_msddelta(logical(roi_I)),'k.')
view(30,30)
grid on
xlabel('MD'), ylabel('var size'), zlabel('mean shape')



% histdat = I_mdiso(logical(roi_I));
% axh_hist = axes('position',[0.6 0.6 0.4 0.4]);
% histogram(axh_hist,histdat)
% 

