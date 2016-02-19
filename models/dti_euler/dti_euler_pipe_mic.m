function fn = dti_euler_pipe_mic(s, o_path, opt)
% function fn = dti_euler_pipe_mic(s, o_path, opt)
%
% s      - input structure
% o_path - output path

if (nargin < 3), opt.present = 1; end
opt = mdm_opt(opt);

% Prepare: mask et c
s = mio_mask_mic(s, o_path, opt);

% Run the analysis
out_fn = fullfile(o_path, 'dti_euler.mat');
mfs_fn = dti_euler_4d_data2fit(s, out_fn, opt);

% Save parameter maps
fn = dti_euler_4d_fit2param(mfs_fn, o_path);

return
    [Itd1,h] = mdm_nii_read(s.nii_fn); 
    xps = s.xps;
    

    Itd1 = squeeze(Itd1(:,:,1,:));
    [nudim.i,nudim.j,td1] = size(Itd1);
    r.i = h.pixdim(2)*(1:nudim.i);
    r.i = r.i - r.i(nudim.i/2+1);
    r.j = h.pixdim(3)*(1:nudim.j);
    r.j = r.j - r.j(nudim.j/2+1);

    Imax = abs(Itd1);
    Imax(nudim.i/2+(0:2),nudim.j/2+(0:2),2:td1) = 0;
    Imax(:,:,1:(td1start-1)) = 0;
    Imax = max(reshape(Imax,numel(Imax),1));                


    Images_S0 = zeros(nudim.i,nudim.j);
    Images_lambda_x = zeros(nudim.i,nudim.j);
    Images_lambda_y = zeros(nudim.i,nudim.j);
    Images_lambda_z = zeros(nudim.i,nudim.j);
    Images_A_alpha = zeros(nudim.i,nudim.j);
    Images_A_beta = zeros(nudim.i,nudim.j);
    Images_A_gamma = zeros(nudim.i,nudim.j);
    Images_Dxx = zeros(nudim.i,nudim.j);
    Images_Dyy = zeros(nudim.i,nudim.j);
    Images_Dzz = zeros(nudim.i,nudim.j);
    Images_lambda1 = zeros(nudim.i,nudim.j);
    Images_lambda2 = zeros(nudim.i,nudim.j);
    Images_lambda3 = zeros(nudim.i,nudim.j);
    Images_v1_x = zeros(nudim.i,nudim.j); Images_v1_y = zeros(nudim.i,nudim.j); Images_v1_z = zeros(nudim.i,nudim.j);
    Images_v2_x = zeros(nudim.i,nudim.j); Images_v2_y = zeros(nudim.i,nudim.j); Images_v2_z = zeros(nudim.i,nudim.j);
    Images_v3_x = zeros(nudim.i,nudim.j); Images_v3_y = zeros(nudim.i,nudim.j); Images_v3_z = zeros(nudim.i,nudim.j);

    bT_xx = xps.bt(:,1);
    bT_yy = xps.bt(:,2);
    bT_zz = xps.bt(:,3);
    bT_xy = xps.bt(:,4)/sqrt(2);
    bT_xz = xps.bt(:,5)/sqrt(2);
    bT_yz = xps.bt(:,6)/sqrt(2);
    bT_trace = xps.b;

    nudim_i = nudim.i;
    nudim_j = nudim.j;

    opt.dti_euler.lsq_opts = optimset('MaxFunEvals',1e4,'Display','off');
    p =  TimedProgressBar( nudim.j, 10, ...
    'Computing. Remaining time: ', ', Completed: ', 'Concluded in ' );

    parfor nj = 1:nudim_j                        
        for ni = 1:nudim_i
            %ni = 10; nj = 10;
            Itd1_pixel = squeeze(Itd1(ni,nj,:));
            ind = td1start:td1;
            signal = abs(Itd1_pixel);
            Iplot = abs(squeeze(Itd1(:,:,opt.mask.b0_ind)));
            Iplot = Iplot/max(max(Iplot));
            if Iplot(ni,nj)>thresh

                m = dti_euler_1d_data2fit(signal, xps, opt, ind);
                                
                Images_S0(ni,nj) = m(1);
                Images_lambda_x(ni,nj) = m(2);
                Images_lambda_y(ni,nj) = m(3);
                Images_lambda_z(ni,nj) = m(4);
                Images_A_alpha(ni,nj) = m(5);
                Images_A_beta(ni,nj) = m(6);
                Images_A_gamma(ni,nj) = m(7);

           end                  
        end
        p.progress; %Counter for progress report
    end
    p.stop;

    Images.dt.s0 = Images_S0;
    Images.dt.lambda.x = Images_lambda_x;
    Images.dt.lambda.y = Images_lambda_y;
    Images.dt.lambda.z = Images_lambda_z;
    Images.dt.euler.alpha = Images_A_alpha;
    Images.dt.euler.beta = Images_A_beta;
    Images.dt.euler.gamma = Images_A_gamma;

if exist([o_path '/dti_euler']) ~= 7
    mkdir(o_path,'dti_euler')
end

Images.r = r;

clear mfs
mfs = Images;
mfs.s = s;
mfs.nii_h = h;

fn =[o_path '/dti_euler/xfs.mat'];
mdm_mfs_save(mfs, s, fn)            

param = {'dt.s0','dt.lambda.x','dt.lambda.y','dt.lambda.z',...
    'dt.euler.alpha','dt.euler.beta','dt.euler.gamma'};
for n = 1:length(param)
    x_fn = [o_path '/dti_euler/' param{n} '.nii'];
    eval(['x = mfs.' param{n} ';'])
    mdm_nii_write(x, x_fn, h);
end
            
