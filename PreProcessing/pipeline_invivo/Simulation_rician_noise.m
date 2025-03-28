%% simulation rician bias correction
clearvars; close all; clc;

addpath('supplementary_functions')
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

%% parameters
image_size = [256 256];
signal_level = [30 100; 30 5;10 2; 10 0];
n_img = 300;
noise_level = 1;

% create images
imgR = ones(image_size(1),image_size(2),n_img);
imgI = zeros(image_size(1),image_size(2),n_img);

% create noise
noiseR = randn(image_size(1),image_size(2),n_img).*noise_level;
noiseI = randn(image_size(1),image_size(2),n_img).*noise_level;

% impart intensity to image
img_int = rand(n_img,1);
[~,sort_order] = sort(img_int);

Nserie = size(signal_level,1);
for ind_serie = 1:Nserie
    disp_int = abs(signal_level(ind_serie,2)-signal_level(ind_serie,1));
    img_intSerie = (img_int.*disp_int)+min(signal_level(ind_serie,:));
    ref_img_int = repmat(permute(img_intSerie,[2 3 1]),image_size(1),image_size(2),1);
    imgRs = imgR.*reshape(img_intSerie,1,1,n_img);
    imgIs = imgI.*reshape(img_intSerie,1,1,n_img);
    
    % create complex image
    img_cmplx = complex(imgRs+imgIs);
    noise_cmplx = complex(noiseR,noiseI);
    
    imgN_cmplx = img_cmplx + noise_cmplx;
    imgN_abs = abs(imgN_cmplx);
    
    %% save to nifti for denoinsing
    imgN_abs = reshape(imgN_abs,[1 2 4 3]);
    my_save_NIFTI(data,[pwd filesep 'data_tmp.nii.gz'])
    
    % compute noise distribution
    NoiseDtmp = imgN_abs-ref_img_int;
    NoiseD(:,ind_serie) = NoiseDtmp(:);
end

% display
edges = [-4*noise_level:(8*noise_level)/1000:4*noise_level]';
x_edges = (edges(1:end-1)+edges(2:end))/2;
figure(1)
for ind_serie = 1:Nserie
    subplot(Nserie,1,ind_serie)
    h=histogram(NoiseD(:,ind_serie),edges);
    Perfect_Gauss = normpdf(x_edges,0,noise_level);
    Perfect_Gauss = Perfect_Gauss./max(Perfect_Gauss);
    Perfect_Gauss = Perfect_Gauss.*max(h.BinCounts);
    hold on
    plot(x_edges,Perfect_Gauss,'r')
    title(['Image serie SNR ' num2str(signal_level(ind_serie,1)) ' to ' num2str(signal_level(ind_serie,2))])
end
% linkaxes


function my_save_NIFTI(data,nii_fn)
% make nifti header
h = mdm_nii_h_empty;
sdim = size(data);
h.pixdim(1+(1:length(sdim))) = sdim;
h.pixdim(2:4) = [1 1 1];
h.xyzt_units = 'SI';
h.dim_info = [1 2 3];
h.sform_code = 0;
mdm_nii_write(data, nii_fn, h, 0);
end



