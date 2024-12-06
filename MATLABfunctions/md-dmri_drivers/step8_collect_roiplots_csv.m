clear all


data_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/processing_dicom/all_RWI_MRI/dicom_100419';
data_dir = dir(fullfile(data_path,'*'));
rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois_20200320';
% rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois_DT';
% rois_path = '/Users/daniel/Dropbox/NMRdata/United/RWI_Matlab/UIH_Wuhan_Brain_Glioma/rois_20200304';

hp_fn = '/Users/daniel/Dropbox/NMRdata/United/histopathology_20200317_DT.csv';

hp_table = readtable(hp_fn);
ind_glioma = hp_table.glioma==1;
data_names = hp_table.case_name(ind_glioma);
grades = hp_table.grade(ind_glioma);

% [~,rois_name,~] = fileparts(rois_path);
% data_names = cell(size(data_dir,1),1);
% for ndata = 1:size(data_dir,1)
%     data_names{ndata} = data_dir(ndata).name;
% end
% data_names = data_names(~contains(data_names,'.'));

roi_prefix = 'lesion';

project_path = fileparts(data_path);
roiplots_path = fullfile(project_path,'roiplots_gliomas');
plot_types = {'dtd_s2000'; 'histograms'; 'technicolor'};
for ntype = 1:numel(plot_types)
    msf_mkdir(fullfile(roiplots_path,plot_types{ntype}))
end

for ndata = 1:numel(data_names)
    data_name = data_names{ndata}; 
    grade = grades(ndata);
    rwi_path = fullfile(data_path,data_name,'rwi');
%     roi_path = fullfile(rois_path,data_name);
        
    roi_dir = dir(fullfile(rwi_path,[roi_prefix '*.pdf']));
    if isempty(roi_dir), continue, end
    
    for npdf = 1:numel(roi_dir)
        in_fn = fullfile(roi_dir(npdf).folder,roi_dir(npdf).name);
        
        for ntype = 1:numel(plot_types)
            plot_type = plot_types{ntype};
            if strfind(roi_dir(npdf).name,plot_type)
                out_fn = fullfile(roiplots_path,plot_type,['grade' num2str(grade) '_' data_name '_' roi_dir(npdf).name]);
                res = copyfile(in_fn,out_fn);
            end
        end
            
    end   
        
    
end
