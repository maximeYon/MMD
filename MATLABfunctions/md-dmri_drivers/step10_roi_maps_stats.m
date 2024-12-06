clear all

data_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom_manmask/all_RWI_MRI/dicom_100419';
% % rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois/DT';
rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois/Yuan20200320';
% rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois/Yuan20200304';
% rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois/Yuan20200415';
hp_fn = '/Users/daniel/Dropbox/NMRdata/United/histopathology_glioma_20200404.xlsx';

% data_path = '/Users/daniel/Dropbox/NMRdata/GE_MSKCC_breast/processing_dicom_manmask/acquisition_date_time';
% rois_path = '/Users/daniel/Dropbox/NMRdata/GE_MSKCC_breast/processing_dicom_manmask/rois_DT';
% hp_fn = '/Users/daniel/Dropbox/NMRdata/GE_MSKCC_breast/MSKCC_histology_IDN20200416_DT.xlsx';

if exist('hp_fn','var')
    opts = detectImportOptions(hp_fn);
    try
        opts.SelectedVariableNames = {'Grade';'CaseID'};
        hp_table = readtable(hp_fn,opts);
        caseids = hp_table.CaseID;
        pathologylabels = cellstr(num2str(hp_table.Grade));
        datalabels = hp_table.CaseID;
    catch
    end
    try
        opts.SelectedVariableNames = {'AcquisitionDateTime';'PathologyLabel';'MolecularLabel';'ScannerStudyID';'Grade'};
        hp_table = readtable(hp_fn,opts);
        caseids = hp_table.AcquisitionDateTime;
        hp_table.PathologyMolecularLabel  = strcat(hp_table.PathologyLabel, '; ', hp_table.MolecularLabel);
        datalabels = hp_table.ScannerStudyID;
        pathologylabels = hp_table.PathologyLabel;
    catch
    end

    isdefined_grades = 1;
else
    isdefined_grades = 0;
end

% Make paths to nii_xps, boostraps, and maps folders
method = 'dtd';

data_dir = dir(fullfile(rois_path,'*'));
[~,rois_name,~] = fileparts(rois_path);
data_names = cell(size(data_dir,1),1);
for ndata = 1:size(data_dir,1)
    data_names{ndata} = data_dir(ndata).name;
end
data_names = data_names(~contains(data_names,'.'));
% data_names = data_names(71);
% data_names = data_names(contains(data_names,{'20191002184711'}));

pmaps_paths = cell(0,0);
for ndata = 1:numel(data_names)
    pmaps_paths{1+numel(pmaps_paths)} = fullfile(data_path,data_names{ndata},'rwi','pmaps');
end

[roi_plots_path,roi_plots_name,~] = fileparts(rois_path);

roi_prefix = 'lesion';
% roi_prefix = 'wm';

% pmap_names = {'dtd_s0'; 'dtd_s1000'; 'dtd_s2000'; 'dti_euler_FA_u_rgb'; 'dti_euler_MD';...
%     'dtd_mdiso';'dtd_msddelta';'dtd_vdiso';'dtd_vsddelta';...
%     'dtd_mdiso_bin1'; 'dtd_mdiso_bin2'; 'dtd_mdiso_bin3';...
%     'dtd_msddelta_bin1'; 'dtd_msddelta_bin2'; 'dtd_msddelta_bin3';...
%     'dtd_mdii_bin1'; 'dtd_mdii_bin2'; 'dtd_mdii_bin3';...
%     'dtd_fractions'};
% roi_plots_path = fullfile(fileparts(roi_plots_path),'roi_plots',[roi_plots_name '_3slice_hist_allmaps']);

% pmap_names = {'dtd_s800'; 'dtd_adc800';...
%     'dtd_msddelta';'dtd_fractions'};
% roi_plots_path = fullfile(fileparts(roi_plots_path),'roi_plots',[roi_plots_name '_3slice_hist_fewmaps_trimmed']);
% 
% pmap_names = {'dtd_s2000'; 'dti_euler_FA_u_rgb'; 'dtd_MD';...
%     'dtd_msddelta';'dtd_fractions'};
% roi_plots_path = fullfile(fileparts(roi_plots_path),'roi_plots',[roi_plots_name '_3slice_hist_fewmaps']);
pmap_names = {'dtd_s2000'; 'dti_euler_FA_u_rgb'; 'dti_euler_MD';...
    'dtd_gamma_MKa';'dtd_gamma_MKi';'dtd_gamma_MKt';...
    'dtd_msddelta';'dtd_vdiso';'dtd_vsddelta';'dtd_fractions';...
    'dtd_mdiso_bin1_gray';'dtd_mdiso_bin2_gray';'dtd_mdiso_bin3_gray';'dtd_mdiso_bin4_gray';'dtd_mdiso_bin5_gray';...
    'dtd_msddelta_bin1_gray';'dtd_msddelta_bin2_gray';'dtd_msddelta_bin3_gray';'dtd_msddelta_bin4_gray';'dtd_msddelta_bin5_gray'};
roi_plots_path = fullfile(fileparts(roi_plots_path),'roi_plots',[roi_plots_name '_3slice_hist_RSNA']);
% 
% pmap_names = {'dtd_s0'; 'dtd_s800'; 'dtd_adc800'};
% roi_plots_path = fullfile(fileparts(roi_plots_path),'roi_plots',[roi_plots_name '_3slice_hist_fewermaps_trimmed']);

% pmap_names = {'dtd_s0'; 'dtd_s1000'; 'dtd_s2000'; ...
%     'dtd_mdiso';'dtd_msddelta';'dtd_vdiso';'dtd_vsddelta';...
%     'dti_euler_FA_u_rgb'};
% roi_plots_path = fullfile(fileparts(roi_plots_path),'roi_plots',[roi_plots_name '_stats']);

% pmap_names = {'dtd_fractions';...
%     'dtd_mdiso_bin1'; 'dtd_mdiso_bin2'; 'dtd_mdiso_bin3';...
%     'dtd_msddelta_bin1'; 'dtd_msddelta_bin2'; 'dtd_msddelta_bin3';...
%     'dtd_mdii_bin1'; 'dtd_mdii_bin2'; 'dtd_mdii_bin3'};
% roi_plots_path = fullfile(fileparts(roi_plots_path),'roi_plots',[roi_plots_name '_bins']);

pmap_ext = '.nii.gz';
msf_mkdir(roi_plots_path)

lw = 1;
fs = 8;
lf_color = [1 1 .999];

Nmaps = numel(pmap_names);

% Define color limits for parameter maps
clim.mdiso = 3.0e-9*[0 1];
clim.adc800 = clim.mdiso/1e-9;
clim.MD = clim.mdiso/1e-9;
clim.MKa = 4*[0 1]; clim.MKi = clim.MKa; clim.MKt = clim.MKa;
clim.msddelta = 1*[0 1];
clim.vdiso = .3*3e-9^2*[0 1];
clim.vsddelta = .15*[0 1];
clim.cvdisosddelta = .1*3e-9*1*[-1 1];
clim.mask_threshold = .02;
Nbins = 20;
fh = figure('Visible', 'off');
% fh = figure('Visible', 'on');
set(fh,'Color','k','InvertHardCopy','off');


table_fn = fullfile(fileparts(hp_fn),[rois_name '_stats.xlsx']);
quantiles_percent = [10 25 50 75 90];
stats_struct.AcquisitionDateTime = data_names;
emptycell = cell(numel(data_names),1);
stats_struct.ROIVolume = emptycell;
try
    stats_struct.ScannerStudyID = hp_table.ScannerStudyID;
catch
end
stats_struct.PathologyLabel = emptycell;
for nmap = 1:Nmaps
    pmap_name = pmap_names{nmap};
    for nquantile = 1:numel(quantiles_percent)
        if strcmp(pmap_name,'dtd_fractions')
            for nbin = 1:3
                column_name = strcat(pmap_name,'_',num2str(nbin),'_',num2str(quantiles_percent(nquantile)));
                stats_struct.(column_name) = emptycell;
            end        
        else   
            if strcmp(pmap_name((end-5):end),'_u_rgb')
                pmap_name = pmap_name(1:(end-6));
            end
            column_name = strcat(pmap_name,'_',num2str(quantiles_percent(nquantile)));
            stats_struct.(column_name) = emptycell;
        end
    end
end

for ndata = 1:numel(data_names)
    
    data_name = data_names{ndata};
    roi_path = fullfile(rois_path,data_name);
    pmaps_path = pmaps_paths{ndata};    
            
    roi_dir = dir(fullfile(roi_path,[roi_prefix '*']));    
    
    if exist(pmaps_path,'dir') ~= 7, disp([data_name ' pmaps_path missing']), continue, end
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
           
    if isdefined_grades    
        ind_table = find(strcmp(caseids,data_name));
%         grade = grades(ind_table);
%         data_name_str = ['grade ' num2str(grade) ' ' data_name];
    else
        data_name_str = [data_name];
    end
    if isempty(ind_table), disp([data_name ' caseid missing']), continue, end
    
    pathologylabel = pathologylabels{ind_table};
    datalabel = datalabels{ind_table};
    stats_struct.PathologyLabel{ndata} = pathologylabel;        
    data_name_str = [pathologylabel '; ' datalabel];

    fig_fn = fullfile(roi_plots_path,[data_name_str '.pdf']);

    opt = mdm_opt();
    opt = mplot_opt(opt);
    opt.roiplot.data_name = data_name;

    try
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


        pmap_name = pmap_names{1};
        pmap_fn = fullfile(pmaps_path,[pmap_name pmap_ext]);

        [I,nii_h] = mdm_nii_read(pmap_fn);
        I = double(I);
        sz = size(I);

        fov.x = nii_h.dim(2)*nii_h.pixdim(2);
        fov.y = nii_h.dim(3)*nii_h.pixdim(3);
        fov.z = nii_h.dim(4)*nii_h.pixdim(4);
        
        papersize = 8.3*[1*Nmaps (fov.y+fov.z)/(fov.y+fov.x)];
        set(fh, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 

        position_hist_0   = [fov.x/(fov.x+fov.y)+.03 .2 fov.y/(fov.x+fov.y)-.18 fov.y/(fov.y+fov.z)-.35];
        position_pmap_x_0 = [fov.x/(fov.x+fov.y) fov.y/(fov.y+fov.z) fov.y/(fov.x+fov.y) fov.z/(fov.y+fov.z)];
        position_pmap_y_0 = [0 fov.y/(fov.y+fov.z) fov.x/(fov.x+fov.y) fov.z/(fov.y+fov.z)];
        position_pmap_z_0 = [0 0 fov.x/(fov.x+fov.y) fov.y/(fov.y+fov.z)];


        stats_struct.ROIVolume{ndata} = sum(roi_I,'all');
        
        [I,nii_h] = mdm_nii_read(fullfile(pmaps_path,['dtd_s0' pmap_ext]));        
        s0_nk = I(:,:,roi_maxz);
        smax = quantile(s0_nk(s0_nk>0),.999,'all');
        mask = I > clim.mask_threshold*smax;
        
        for nmap = 1:Nmaps
            pmap_name = pmap_names{nmap};
            pmap_fn = fullfile(pmaps_path,[pmap_name pmap_ext]);
            if exist(pmap_fn,'file') ~= 2, disp([data_name ' ' pmap_name ' missing']), continue, end

            npmap = nmap;

            position_hist     = position_hist_0.*[1/Nmaps 1 1/Nmaps 1] + [(npmap-1)/Nmaps 0 0 0];
            position_pmap_x   = position_pmap_x_0.*[1/Nmaps 1 1/Nmaps 1] + [(npmap-1)/Nmaps 0 0 0];
            position_pmap_y   = position_pmap_y_0.*[1/Nmaps 1 1/Nmaps 1] + [(npmap-1)/Nmaps 0 0 0];
            position_pmap_z   = position_pmap_z_0.*[1/Nmaps 1 1/Nmaps 1] + [(npmap-1)/Nmaps 0 0 0];


            [I,nii_h] = mdm_nii_read(pmap_fn);
            I = double(I);
            sz = size(I);            

            axh_hist = axes('position',position_hist,'Color','k','xcolor',lf_color);
            hold(axh_hist,'on')


            if sz(1) == 3
                roi_col = [1 1 .999];
                clim_temp = [0 1];
                edges = linspace(clim_temp(1),clim_temp(2),Nbins);
                sz = sz(2:4);
                I_r = reshape(I(1,:,:,:),sz)/255.*mask;
                I_g = reshape(I(2,:,:,:),sz)/255.*mask;
                I_b = reshape(I(3,:,:,:),sz)/255.*mask; 

                im2d_x = cat(3,squeeze(I_r(roi_maxx,:,:))',squeeze(I_g(roi_maxx,:,:))',squeeze(I_b(roi_maxx,:,:))');
                im2d_y = cat(3,squeeze(I_r(:,roi_maxy,:))',squeeze(I_g(:,roi_maxy,:))',squeeze(I_b(:,roi_maxy,:))');
                im2d_z = cat(3,squeeze(I_r(:,:,roi_maxz))',squeeze(I_g(:,:,roi_maxz))',squeeze(I_b(:,:,roi_maxz))');

                histdat_r = I_r;        
                histdat_g = I_g;        
                histdat_b = I_b;        
                if strcmp(pmap_name((end-8):end),'fractions')
                    roi_col = [0 0 0];
                    fractions_norm = sum(cat(4,histdat_r,histdat_g,histdat_b),4);
                    histdat_r = histdat_r./fractions_norm;        
                    histdat_g = histdat_g./fractions_norm;        
                    histdat_b = histdat_b./fractions_norm;  
                    bin_col = {'r';'g';'b'};
                    for nquantile = 1:numel(quantiles_percent)
                        for nbin = 1:3
                            eval(['histdat = histdat_' bin_col{nbin} ';'])
                            quantiles_val = quantile(histdat(roi_I),quantiles_percent(nquantile)/100);
                            column_name = strcat(pmap_name,'_',num2str(nbin),'_',num2str(quantiles_percent(nquantile)));
                            stats_struct.(column_name){ndata} = quantiles_val;
                        end
                    end
                end
                [counts_r,~] = histcounts(histdat_r(roi_I),edges);
                [counts_g,~] = histcounts(histdat_g(roi_I),edges);
                [counts_b,~] = histcounts(histdat_b(roi_I),edges);
                hh_r = histogram(axh_hist,'BinEdges', edges, 'BinCounts', counts_r);
                hh_g = histogram(axh_hist,'BinEdges', edges, 'BinCounts', counts_g);
                hh_b = histogram(axh_hist,'BinEdges', edges, 'BinCounts', counts_b);
                ylim = max([counts_r(2:(end-1)) counts_g(2:(end-1)) counts_b(2:(end-1))])*[-.1 1.1];
                set(hh_r,'DisplayStyle','stairs','EdgeColor','r')
                set(hh_g,'DisplayStyle','stairs','EdgeColor','g')
                set(hh_b,'DisplayStyle','stairs','EdgeColor','b')

                if strcmp(pmap_name((end-5):end),'_u_rgb')
                    delete([hh_r; hh_g; hh_b])
                    %hold(axh_hist,'off')
                    pmap_name = pmap_name(1:(end-6));
                    pmap_fn = fullfile(pmaps_path,[pmap_name pmap_ext]);

                    [histdat,~] = mdm_nii_read(pmap_fn);
                    histdat = double(histdat);
                    [counts,~] = histcounts(histdat(roi_I),edges);
                    hh = histogram(axh_hist,'BinEdges', edges, 'BinCounts', counts);
                    ylim = max(counts(2:(end-1)))*[-.1 1.1];
                    set(hh,'DisplayStyle','stairs','EdgeColor',lf_color)
                    for nquantile = 1:numel(quantiles_percent)
                        quantiles_val = quantile(histdat(roi_I),quantiles_percent(nquantile)/100);
                        column_name = strcat(pmap_name,'_',num2str(quantiles_percent(nquantile)));
                        stats_struct.(column_name){ndata} = quantiles_val;
                    end
                end
            else
                roi_col = [1 0 0];
                p_name = pmap_name((max(find(pmap_name=='_'))+1):end);
                if isfield(clim,p_name)
                    clim_temp = clim.(p_name);
                else
                    I_maxz = I(:,:,roi_maxz);
                    clim_temp = quantile(I_maxz(I_maxz>0),.999,'all')*[0 1];
                end
                edges = linspace(clim_temp(1),clim_temp(2),Nbins);

                I_clamp = mio_min_max_cut(I, clim_temp(1), clim_temp(2)).*mask;
                I_r = I_clamp;
                I_g = I_clamp;
                I_b = I_clamp;   

                histdat = I;
                [counts,~] = histcounts(histdat(roi_I),edges);
                hh = histogram(axh_hist,'BinEdges', edges, 'BinCounts', counts);
                ylim = max(counts(2:(end-1)))*[-.1 1.1];
                set(hh,'DisplayStyle','stairs','EdgeColor',lf_color)
                for nquantile = 1:numel(quantiles_percent)
                    quantiles_val = quantile(histdat(roi_I),quantiles_percent(nquantile)/100);
                    column_name = strcat(pmap_name,'_',num2str(quantiles_percent(nquantile)));
                    stats_struct.(column_name){ndata} = quantiles_val;
                end
            end

            im2d_x = cat(3,squeeze(I_r(roi_maxx,:,:))',squeeze(I_g(roi_maxx,:,:))',squeeze(I_b(roi_maxx,:,:))');
            im2d_y = cat(3,squeeze(I_r(:,roi_maxy,:))',squeeze(I_g(:,roi_maxy,:))',squeeze(I_b(:,roi_maxy,:))');
            im2d_z = cat(3,squeeze(I_r(:,:,roi_maxz))',squeeze(I_g(:,:,roi_maxz))',squeeze(I_b(:,:,roi_maxz))');

            axh_pmap_x = axes('position',position_pmap_x);
            axh_pmap_y = axes('position',position_pmap_y);
            axh_pmap_z = axes('position',position_pmap_z);

            imagesc(axh_pmap_x,.99*im2d_x/clim_temp(2)) % Factor .99 to avoid black pixels in pdf
            imagesc(axh_pmap_y,.99*im2d_y/clim_temp(2))
            imagesc(axh_pmap_z,.99*im2d_z/clim_temp(2))

            set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'YDir','normal')

            plot_roi(roi2d_x',roi_col, .5*lw, 1, axh_pmap_x);
            plot_roi(roi2d_y',roi_col, .5*lw, 1, axh_pmap_y);
            plot_roi(roi2d_z',roi_col, .5*lw, 1, axh_pmap_z);

            set([axh_pmap_x; axh_pmap_y; axh_pmap_z],'CLim',clim_temp)
            axis([axh_pmap_x; axh_pmap_y; axh_pmap_z],'off')


            set(axh_hist,'XLim',clim_temp,'YLim',ylim)
            set(axh_hist,'Box','off','TickDir','out','TickLength',.02*[1 1],'LineWidth',lw,'FontSize',fs,'YTick',[])
            xlabel(axh_hist,pmap_name,'FontSize',fs,'Interpreter','none')

            if nmap == 1
                hold(axh_pmap_x,'on')
                hold(axh_pmap_y,'on')
                hold(axh_pmap_z,'on')
                phxz = plot(axh_pmap_z,roi_maxx*[1 1],[-1 sz(1)+1],'-');
                phyz = plot(axh_pmap_z,[-1 sz(2)+1],roi_maxy*[1 1],'-');
                phzy = plot(axh_pmap_y,[-1 sz(2)+1],roi_maxz*[1 1],'-');
                phxy = plot(axh_pmap_y,roi_maxx*[1 1],[-1 sz(3)+1],'-');
                phzx = plot(axh_pmap_x,[-1 sz(1)+1],roi_maxz*[1 1],'-');
                phyx = plot(axh_pmap_x,roi_maxy*[1 1],[-1 sz(3)+1],'-');
                set([phzy; phzx; phxy; phxz; phyz; phyx],'Color',[1 1 .999])

                th_x = text(axh_pmap_x,.5,sz(3)+.5,['x = ' num2str(roi_maxx)]);
                th_y = text(axh_pmap_y,.5,sz(3)+.5,['y = ' num2str(roi_maxy)]);
                th_z = text(axh_pmap_z,.5,sz(2)+.5,['z = ' num2str(roi_maxz)]);
                set([th_x; th_y; th_z],'Color',[1 1 .999],'VerticalAlignment','top')
                th_data_z = text(axh_pmap_z,.5,.5,data_name_str);
                set(th_data_z,'Color',[1 1 .999],'VerticalAlignment','bottom')
                set([th_z; th_y; th_x; th_data_z],'Color',[1 1 .999],'FontSize',fs)
            end

        end

        print(fh,fig_fn,'-loose','-dpdf')     
        clf(fh)
    catch me
        disp([data_name])
        disp( getReport( me, 'extended', 'hyperlinks', 'on' ) )
    end       
end

stats_table = struct2table(stats_struct);
writetable(stats_table,table_fn)


