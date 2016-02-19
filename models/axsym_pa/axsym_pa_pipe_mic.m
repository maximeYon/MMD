function fn = axsym_pa_pipe_mic(s, o_path, opt)
% function fn = dti_nls_pipe_example(s, o_path, opt)
%
% s      - input structure
% o_path - output path

if (nargin < 3), opt.present = 1; end
opt = mdm_opt(opt);

td1start = opt.td1start;
thresh = opt.mask.thresh;

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
    Images_S0_gamma = zeros(nudim.i,nudim.j);
    Images_MD_gamma = zeros(nudim.i,nudim.j);
    Images_mu2i = zeros(nudim.i,nudim.j);
    Images_mu2a_gamma = zeros(nudim.i,nudim.j);
    Images_Ki = zeros(nudim.i,nudim.j);
    Images_Ka_gamma = zeros(nudim.i,nudim.j);
    Images_uFA_gamma = zeros(nudim.i,nudim.j);
    Images_S0_erf = zeros(nudim.i,nudim.j);
    Images_MD_erf = zeros(nudim.i,nudim.j);
    Images_DDelta = zeros(nudim.i,nudim.j);
    Images_mu2a_erf = zeros(nudim.i,nudim.j);
    Images_Ka_erf = zeros(nudim.i,nudim.j);
    Images_uFA_erf = zeros(nudim.i,nudim.j);

    Images_S_PA = zeros(nudim.i,nudim.j,xps.b_nind*xps.bd_nind);
    Images_S_PA_gamma = zeros(nudim.i,nudim.j,xps.b_nind*xps.bd_nind);
    Images_S_PA_erf_oblate = zeros(nudim.i,nudim.j,xps.b_nind*xps.bd_nind);
    Images_S_PA_erf_prolate = zeros(nudim.i,nudim.j,xps.b_nind*xps.bd_nind);


    bT_xx = xps.bt(:,1);
    bT_yy = xps.bt(:,2);
    bT_zz = xps.bt(:,3);
    bT_xy = xps.bt(:,4)/sqrt(2);
    bT_xz = xps.bt(:,5)/sqrt(2);
    bT_yz = xps.bt(:,6)/sqrt(2);
    bT_trace = xps.b;

    nudim_i = nudim.i;
    nudim_j = nudim.j;

    options = optimset('MaxFunEvals',1e4,'Display','off');
    p =  TimedProgressBar( nudim.j, 10, ...
    'Computing. Remaining time: ', ', Completed: ', 'Concluded in ' );

    parfor nj = 1:nudim_j                        
        for ni = 1:nudim_i
            Itd1_pixel = squeeze(Itd1(ni,nj,:));
            indx_fit = td1start:td1;
            PA_weight = zeros(xps.n,1); PA_weight(indx_fit) = 1;
            PAweight_array = reshape(PA_weight,xps.b_nind*xps.bd_nind,xps.br_nind);
            S_array = reshape(Itd1_pixel,xps.b_nind*xps.bd_nind,xps.br_nind);
            S_PA = abs(sum(S_array.*PAweight_array,2)./sum(PAweight_array,2));

            %figure(1), clf, plot(1:NbTtrace*NbTDelta,abs(S_PA)), return

            %Yin = squeeze(sum(sum(abs(Itd1(:,:,indx_fit)),1),2));
            Yin = abs(Itd1_pixel(indx_fit));
            Iplot = abs(squeeze(Itd1(:,:,opt.mask.b0_ind)));
            Iplot = Iplot/max(max(Iplot));
            if Iplot(ni,nj)>thresh
                Images_S_PA(ni,nj,:) = S_PA;

                Xin_xx = bT_xx(indx_fit);
                Xin_xy = bT_xy(indx_fit);
                Xin_xz = bT_xz(indx_fit);
                Xin_yy = bT_yy(indx_fit);
                Xin_yz = bT_yz(indx_fit);
                Xin_zz = bT_zz(indx_fit);

                Pin = [max(Yin)   1.5e-9 1e-9 .5e-9 2*pi 2*pi 2*pi]; Funam = 'fexpDTImat';
                LB = [0 1e-11*[1 1 1] 0 0 0]; UB = [2*max(Yin) 10e-9*[1 1 1] 4*pi 4*pi 4*pi];

                %figure(1), clf, plot(indx_fit,Yin,'o'), return
                Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin_xx); Pnorm = abs(Pin); 
                [Pout,resnorm,residual,exitflag,output] = lsqcurvefit(Funam,Pin./Pnorm,Xin_xx/Xnorm,...
                    Yin/Ynorm,LB./Pnorm,UB./Pnorm,options,...
                    Xin_xy/Xnorm,Xin_xz/Xnorm,Xin_yy/Xnorm,Xin_yz/Xnorm,Xin_zz/Xnorm,...
                    Pnorm,Xnorm,Ynorm); Pout = Pout.*Pnorm;
                Ycalc = feval(Funam,Pout,...
                    Xin_xx,Xin_xy,Xin_xz,Xin_yy,Xin_yz,Xin_zz,...
                    ones(size(Pin)),1,1); error = Yin - Ycalc;
                %return

%                             if PlotInterm
%                                 figure(2), clf
%                                 subplot(2,2,3)
%                                 plot(indx_fit,Yin,'-o',indx_fit,Ycalc,'-x'), title(['ni=' num2str(ni) ' nj=' num2str(nj)]);
%                                 %return
%                             end

                S0 = Pout(1);
                lambda_x = Pout(2);
                lambda_y = Pout(3);
                lambda_z = Pout(4);
                A_alpha = angle(exp(i*Pout(5)));
                A_beta = angle(exp(i*Pout(6)));
                A_gamma = angle(exp(i*Pout(7)));

                R_gamma = [
                    cos(A_gamma) sin(A_gamma) 0
                    -sin(A_gamma) cos(A_gamma) 0
                    0 0 1];

                R_beta = [
                    cos(A_beta) 0 -sin(A_beta)
                    0 1 0
                    sin(A_beta) 0 cos(A_beta)];

                R_alpha = [
                    cos(A_alpha) sin(A_alpha) 0
                    -sin(A_alpha) cos(A_alpha) 0
                    0 0 1];

                R_mat = R_gamma*R_beta*R_alpha;
                R_inv = inv(R_mat);

                DT_PAS = [
                    lambda_x 0 0
                    0 lambda_y 0
                    0 0 lambda_z];

                DT_LF = R_mat*DT_PAS*R_inv;

                lambdas = [lambda_x lambda_y lambda_z];
                [lambdas,indx] = sort(lambdas,2,'descend');
                lambda1 = lambdas(1);
                lambda2 = lambdas(2);
                lambda3 = lambdas(3);

                vx_x = R_mat(1,1); vx_y = R_mat(2,1); vx_z = R_mat(3,1);
                vy_x = R_mat(1,2); vy_y = R_mat(2,2); vy_z = R_mat(3,2);
                vz_x = R_mat(1,3); vz_y = R_mat(2,3); vz_z = R_mat(3,3);

                v1_x = 0; v1_y = 0; v1_z = 0;
                v2_x = 0; v2_y = 0; v2_z = 0;
                v3_x = 0; v3_y = 0; v3_z = 0;
                if lambda1 == lambda_x            
                    if lambda2 == lambda_y
                        v1_x = vx_x; v1_y = vx_y; v1_z = vx_z;
                        v2_x = vy_x; v2_y = vy_y; v2_z = vy_z;
                        v3_x = vz_x; v3_y = vz_y; v3_z = vz_z;
                    elseif lambda2 == lambda_z
                        v1_x = vx_x; v1_y = vx_y; v1_z = vx_z;
                        v2_x = vz_x; v2_y = vz_y; v2_z = vz_z;
                        v3_x = vy_x; v3_y = vy_y; v3_z = vy_z;
                    end
                elseif lambda1 == lambda_y            
                    if lambda2 == lambda_x
                        v1_x = vy_x; v1_y = vy_y; v1_z = vy_z;
                        v2_x = vx_x; v2_y = vx_y; v2_z = vx_z;
                        v3_x = vz_x; v3_y = vz_y; v3_z = vz_z;
                    elseif lambda2 == lambda_z
                        v1_x = vy_x; v1_y = vy_y; v1_z = vy_z;
                        v2_x = vz_x; v2_y = vz_y; v2_z = vz_z;
                        v3_x = vx_x; v3_y = vx_y; v3_z = vx_z;
                    end
                elseif lambda1 == lambda_z            
                    if lambda2 == lambda_x
                        v1_x = vz_x; v1_y = vz_y; v1_z = vz_z;
                        v2_x = vx_x; v2_y = vx_y; v2_z = vx_z;
                        v3_x = vy_x; v3_y = vy_y; v3_z = vy_z;
                    elseif lambda2 == lambda_y
                        v1_x = vz_x; v1_y = vz_y; v1_z = vz_z;
                        v2_x = vy_x; v2_y = vy_y; v2_z = vy_z;
                        v3_x = vx_x; v3_y = vx_y; v3_z = vx_z;
                    end
                end

                Images_S0(ni,nj) = S0;
                Images_lambda_x(ni,nj) = lambda_x;
                Images_lambda_y(ni,nj) = lambda_y;
                Images_lambda_z(ni,nj) = lambda_z;
                Images_A_alpha(ni,nj) = A_alpha;
                Images_A_beta(ni,nj) = A_beta;
                Images_A_gamma(ni,nj) = A_gamma;
                Images_Dxx(ni,nj) = DT_LF(1,1);
                Images_Dyy(ni,nj) = DT_LF(2,2);
                Images_Dzz(ni,nj) = DT_LF(3,3);
                Images_lambda1(ni,nj) = lambda1;
                Images_lambda2(ni,nj) = lambda2;
                Images_lambda3(ni,nj) = lambda3;
                Images_v1_x(ni,nj) = v1_x; Images_v1_y(ni,nj) = v1_y; Images_v1_z(ni,nj) = v1_z;
                Images_v2_x(ni,nj) = v2_x; Images_v2_y(ni,nj) = v2_y; Images_v2_z(ni,nj) = v2_z;
                Images_v3_x(ni,nj) = v3_x; Images_v3_y(ni,nj) = v3_y; Images_v3_z(ni,nj) = v3_z;

                MD = (lambda_x+lambda_y+lambda_z)/3;                            
                mu2_FA = 4/45*((lambda_z-lambda_x)^2+(lambda_y-lambda_z)*(lambda_y-lambda_x));
                FA = sqrt(1/2)*sqrt((lambda_z-lambda_y).^2+(lambda_z-lambda_x).^2+(lambda_y-lambda_x).^2)...
                ./sqrt(lambda_z.^2+lambda_y.^2+lambda_x.^2);

                indx_fit = (1:xps.b_nind*xps.bd_nind);
                Xin = bT_trace(indx_fit);

                ThreshFit = opt.mask.thresh;
                Xthresh = -log(ThreshFit)/MD;
                fitpoints = find(exp(-Xin*MD)>ThreshFit); fitpoints = 1:length(Xin);

                Xin  = xps.b(1:xps.b_nind*xps.bd_nind);
                Xin2  = xps.b_delta(1:xps.b_nind*xps.bd_nind);
                Yin = abs(S_PA);
                Weight = .5-.5*erf(5*(Xin-Xthresh)/Xthresh);
                %figure(1), clf, plot(Xin,Weight,'o'), return
                %figure(1), clf, semilogy(Xin,Yin,'o'), return

                Pin = [S0 1.5*MD (.1*MD)^2 (.4*MD)^2]; Funam = 'fexpgamma_AxSym';
                LB = [.95*S0 0.5*MD (.001*MD)^2 .5*mu2_FA ]; UB = [1.5*S0 2*MD MD^2 MD^2];

                Pout = Pin; Ynorm = mean(Yin)./(Weight+1e-10 ); Xnorm = mean(Xin); Pnorm = Pin; 
                Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin./Ynorm,...
                    LB./Pnorm,UB./Pnorm,options,Pnorm,Xnorm,Ynorm,Xin2);
                Yout = feval(Funam,Pout,Xin,ones(size(Pout)),1,1,Xin2); error = Yin - Yout;
                chisq_gamma = sum(error.^2);
                %figure(1), clf, semilogy(Xin,Yin,'o',Xin,Yout,'-'), return

                S0_gamma = Pout(1);
                MD_gamma = Pout(2);
                mu2i = Pout(3); 
                mu2a_gamma = Pout(4);

                Ycalc_MD = S0_gamma*exp(-MD_gamma*Xin);

                Images_S_PA_gamma(ni,nj,:) = Yout;

                uFA_gamma = sqrt(3/2)*sqrt(1/(1+2/5*MD_gamma^2/(mu2a_gamma)));

                Ki = mu2i/MD_gamma^2;
                Ka_gamma = mu2a_gamma/MD_gamma^2;

%%
                Pin = [S0_gamma MD_gamma sqrt(5/4*mu2a_gamma/MD_gamma^2)]; Funam = 'fDiffVACSY2';
                UB = [max(S0_gamma)*1.1 MD_gamma*2 1];
                LB = [max(S0_gamma)*0.9 MD_gamma*.5 0];

                Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = abs(Pin); 
                Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,LB./Pnorm,UB./Pnorm,options,Xin2,Pnorm,Xnorm,Ynorm);
                Yout = feval(Funam,Pout,Xin,Xin2); error = Yin - Yout;
                chisq_prolate = sum(error.^2);
                Pout_prolate = Pout; Yout_prolate = Yout;

                Pin = [S0_gamma MD_gamma -sqrt(5/4*mu2a_gamma/MD_gamma^2)]; Funam = 'fDiffVACSY2';
                UB = [max(S0_gamma)*1.1 MD_gamma*2 0];
                LB = [max(S0_gamma)*0.9 MD_gamma*.5 -.5];

                Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = abs(Pin); 
                Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,LB./Pnorm,UB./Pnorm,options,Xin2,Pnorm,Xnorm,Ynorm);
                Yout = feval(Funam,Pout,Xin,Xin2); error = Yin - Yout;
                chisq_oblate = sum(error.^2);
                Pout_oblate = Pout; Yout_oblate = Yout;

                semilogy(reshape(Xin,[xps.b_nind xps.bd_nind]),reshape(Yout_prolate,[xps.b_nind xps.bd_nind]),'r-')
                semilogy(reshape(Xin,[xps.b_nind xps.bd_nind]),reshape(Yout_oblate,[xps.b_nind xps.bd_nind]),'b-')
                %figure(1), clf, semilogy(Xin,Yin,'o',Xin,Yout,'-'), return

                Images_S_PA_erf_prolate(ni,nj,:) = Yout_prolate;
                Images_S_PA_erf_oblate(ni,nj,:) = Yout_oblate;

                if chisq_oblate < chisq_prolate
                    Pout = Pout_oblate;
                else
                    Pout = Pout_prolate;
                end
                mu2a_oblate = 4/5*(Pout_oblate(3).*Pout_oblate(2)).^2;
                mu2a_prolate = 4/5*(Pout_prolate(3).*Pout_prolate(2)).^2;

                
                S0_erf = Pout(1);
                MD_erf = Pout(2);
                DDelta = Pout(3);
                mu2a_erf = 4/5*(DDelta.*MD_erf).^2;
                Ka_erf = mu2a_erf/MD_erf^2;
                uFA_erf = sqrt(3/2)*sqrt(1/(1+2/5*MD_erf^2/(mu2a_erf)));
                             
                Images_S0_gamma(ni,nj) = S0_gamma;
                Images_MD_gamma(ni,nj) = MD_gamma;
                Images_mu2i(ni,nj) = mu2i;
                Images_mu2a_gamma(ni,nj) = mu2a_gamma;
                Images_Ki(ni,nj) = Ki;
                Images_Ka_gamma(ni,nj) = Ka_gamma;
                Images_uFA_gamma(ni,nj) = uFA_gamma;
                Images_S0_erf(ni,nj) = S0_erf;
                Images_MD_erf(ni,nj) = MD_erf;
                Images_DDelta(ni,nj) = DDelta;
                Images_mu2a_erf(ni,nj) = mu2a_erf;
                Images_Ka_erf(ni,nj) = Ka_erf;
                Images_uFA_erf(ni,nj) = uFA_erf;

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
    Images.dt.xx = Images_Dxx;
    Images.dt.yy = Images_Dyy;
    Images.dt.zz = Images_Dzz;
    Images.dt.lambda1 = Images_lambda1;
    Images.dt.lambda2 = Images_lambda2;
    Images.dt.lambda3 = Images_lambda3;
    Images.dt.v1.x = Images_v1_x; Images.dt.v1.y = Images_v1_y; Images.dt.v1.z = Images_v1_z;
    Images.dt.v2.x = Images_v2_x; Images.dt.v2.y = Images_v2_y; Images.dt.v2.z = Images_v2_z;
    Images.dt.v3.x = Images_v3_x; Images.dt.v3.y = Images_v3_y; Images.dt.v3.z = Images_v3_z;
    Images.gamma.s0 = Images_S0_gamma;
    Images.gamma.md = Images_MD_gamma;
    Images.gamma.mu2i = Images_mu2i;
    Images.gamma.mu2a = Images_mu2a_gamma;
    Images.gamma.ki = Images_Ki;
    Images.gamma.ka = Images_Ka_gamma;
    Images.gamma.ufa = Images_uFA_gamma;
    Images.erf.s0 = Images_S0_erf;
    Images.erf.md = Images_MD_erf;
    Images.erf.delta = Images_DDelta;
    Images.erf.mu2a = Images_mu2a_erf;
    Images.erf.ka = Images_Ka_erf;
    Images.erf.ufa = Images_uFA_erf;

    Images.s_pa = Images_S_PA;
    Images.gamma.s_pa = Images_S_PA_gamma;
    Images.erf.s_pa.oblate = Images_S_PA_erf_oblate;
    Images.erf.s_pa.prolate = Images_S_PA_erf_prolate;

    Images.dt.md = 1/3*(Images.dt.lambda1+Images.dt.lambda2+Images.dt.lambda3);

    Images.dt.fa = sqrt(1/2)*sqrt((Images.dt.lambda1-Images.dt.lambda2).^2+(Images.dt.lambda1-Images.dt.lambda3).^2+(Images.dt.lambda2-Images.dt.lambda3).^2)...
        ./sqrt(Images.dt.lambda1.^2+Images.dt.lambda2.^2+Images.dt.lambda3.^2);
    Images.dt.fa(find(isnan(Images.dt.fa)==1)) = 0;

    Images.dt.cl = (Images.dt.lambda1 - Images.dt.lambda2)./Images.dt.lambda1;
    Images.dt.cl(find(isnan(Images.dt.cl)==1)) = 0;
    Images.dt.cp = (Images.dt.lambda2 - Images.dt.lambda3)./Images.dt.lambda1;
    Images.dt.cp(find(isnan(Images.dt.cp)==1)) = 0;

if exist([o_path '/axsym']) ~= 7
    mkdir(o_path,'axsym')
end

Images.r = r;

clear mfs
mfs = Images;
mfs.s = s;
mfs.nii_h = h;

fn =[o_path '/axsym/xfs.mat'];
mdm_mfs_save(mfs, s, fn)            

param = {'dt.s0','dt.md','dt.fa','dt.cp','dt.cl',...
    'gamma.s0','gamma.md','gamma.ufa','gamma.mu2i','gamma.mu2a','gamma.ki','gamma.ka',...
    'erf.s0','erf.md','erf.ufa','erf.mu2a','erf.ka','erf.delta'};
for n = 1:length(param)
    x_fn = [o_path '/axsym/' param{n} '.nii'];
    eval(['x = mfs.' param{n} ';'])
    mdm_nii_write(x, x_fn, h);
end
            
