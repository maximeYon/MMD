function out_path = mdm_bruker_dt_moacrare2drs_recon(data_path, out_path, rps, opt)
% function out_path = mdm_bruker_dt_moacrare2drs_recon(data_path, out_path, rps, opt)
%
% Recon images acquired with Bruker pulse sequence DT_moacrare2drs.
%
% Works in progress

if (nargin < 5), opt.present = 1; end

opt = mdm_opt(opt);


% Book-keeping

msf_mkdir(fullfile(out_path)); 
nii_temp_fn = fullfile(out_path, ['data_temp' opt.nii_ext]);
xps_temp_fn = fullfile(out_path, 'data_temp_xps.mat');

% convert bruker acquistion parameters to xps
mdm_bruker_acqus2mat(data_path);
mdm_bruker_dt_moacrare2drs_acqus2xps(data_path, xps_temp_fn);
mdm_bruker_dt_moacrare2drs_ser2nii(data_path, nii_temp_fn, rps);

s.nii_fn = nii_temp_fn;
load(xps_temp_fn);
s.xps = xps;

Nslices = numel(unique(xps.fq1));
NRDenc = xps.n/Nslices;
[~,~,slice_ind] = unique(xps.fq1, 'rows');
for nslice = 1:Nslices
    ind = nslice == slice_ind;
    slice_path = fullfile(fileparts(s.nii_fn),'slices',num2str(nslice));
    s_new.nii_fn = fullfile(slice_path,'nii_xps',['data' opt.nii_ext]);
    
    s_new.xps = mdm_xps_subsample(s.xps, ind);
    mdm_nii_subsample(s.nii_fn, ind, s_new.nii_fn);
    xps = s_new.xps;
    mdm_xps_save(xps, mdm_xps_fn_from_nii_fn(s_new.nii_fn));

%     suffix = ['slice' num2str(nslice)];
%     s_new = mdm_s_subsample(s, ind, fileparts(s.nii_fn), opt, suffix);
end                 

delete(nii_temp_fn);
delete(xps_temp_fn);

delete(fullfile(data_path,'NMRacqus.mat'))
delete(fullfile(data_path,'NMRacqu2s.mat'))
