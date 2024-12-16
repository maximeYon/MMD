function out_path = mdm_bruker_dt_rare2d_recon(data_path, out_path, rps, opt)
% function out_path = mdm_bruker_dt_rare2d_recon(data_path, out_path, rps, opt)
%
% Recon images acquired with Bruker pulse sequence DT_axderare2d.
%
% Works in progress

if (nargin < 5), opt.present = 1; end

opt = mdm_opt(opt);

if isdir(data_path)
    msf_mkdir(fullfile(out_path)); 

    nii_temp_fn = fullfile(out_path, ['data_temp' opt.nii_ext]);
    xps_temp_fn = fullfile(out_path, 'data_temp_xps.mat');

    nii_fn = fullfile(out_path, ['data' opt.nii_ext]);
    xps_fn = fullfile(out_path, 'data_xps.mat');
end

try
    % convert bruker acquistion parameters to xps
    mdm_bruker_acqus2mat(data_path);
    mdm_bruker_dt_axderare2d_acqus2xps(data_path, xps_temp_fn, rps);
    mdm_bruker_dt_rare2d_ser2nii(data_path, nii_temp_fn, rps);

    s.nii_fn = nii_temp_fn;
    load(xps_temp_fn);
    s.xps = xps;

    if isfield(rps,'ind_start')
        ind_start = rps.ind_start;
    else
        ind_start = 2; % Exclude first image
    end

    if ind_start > 1
        mdm_s_subsample(s, (1:s.xps.n) >= ind_start);

        nii_sub_fn = fullfile(out_path, ['data_temp_sub' opt.nii_ext]);
        xps_sub_fn = fullfile(out_path, 'data_temp_sub_xps.mat');

        movefile(nii_sub_fn,nii_fn);
        movefile(xps_sub_fn,xps_fn);

        delete(nii_temp_fn);
        delete(xps_temp_fn);
    else
        movefile(nii_temp_fn,nii_fn);
        movefile(xps_temp_fn,xps_fn);
    end


    delete(fullfile(data_path,'NMRacqus.mat'))
    delete(fullfile(data_path,'NMRacqu2s.mat'))
catch
     disp(['Not processed: ' data_path])
end
