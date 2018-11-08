function color_gamma_plot(S, xps, axh, axh2)
% function color_gamma_plot(S, xps, axh, axh2)

if (nargin < 4), axh2 = []; end

opt = mdm_opt();
opt = dti_euler_opt(opt);
opt = color_gamma_opt(opt);
opt = mplot_opt(opt);
xps = color_gamma_check_xps(xps);

% % Customize options
% opt.dtd.dmin = .2/max(xps.b);

% Fit and predict signal
m     = dti_euler_1d_data2fit(S, xps, opt); 
S_fit = dti_euler_1d_fit2data(m, xps);
m_dti_euler = m;

S0 = m(1);
MD = mean(m(2:4));

s0          = m(1);
lambdax     = m(2);
lambday     = m(3);
lambdaz     = m(4);
euler_alpha = angle(exp(1i*m(5)));
euler_beta  = angle(exp(1i*m(6)));
euler_gamma = angle(exp(1i*m(7)));

[rotmat,rotmatinv] = tm_euler_angles2rotmat(euler_alpha,euler_beta,euler_gamma);

dt_lambda = diag([lambdax, lambday, lambdaz]);

dt3x3 = rotmat*dt_lambda*rotmatinv;
dt = tm_3x3_to_tpars(dt3x3);


% Clear, plot, and configure 
cla(axh); hold(axh, 'off');
cla(axh2); hold(axh2, 'off');


unit_to_SI = [max(S+eps) 1e-9 (1e-9)^2];

m_guess   = [S0 MD .01*MD^2];
m_lb      = [0.99*S0 1e-11 -1e-19];
m_ub      = [1.01*S0 3e-9 5e-18];
                
t_guess   = m_guess./unit_to_SI;
t_lb      = m_lb./unit_to_SI;
t_ub      = m_ub./unit_to_SI;

    function m = t2m(t) % convert local params to outside format
        
        % define model parameters
        s0    = t(1);         % Signal at no diff encoding
        d     = t(2);         % 1st moment
        mu2   = t(3);         % 2nd moment
        
        m = [s0 d mu2] .* unit_to_SI;
    end

    function s = my_1d_fit2data(t, b)
        m = t2m(t);
        
        s0    = m(1);         % Signal at no diff encoding
        d     = m(2);         % 1st moment
        mu2   = m(3);         % 2nd moment
        s = s0.*(1 + b.*mu2./d).^(-1*(d.^2./mu2));
        s = real(s).*weight;
    end
%options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','Display','off');
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display','off');

Smax = max(S(:));
bmax = max(xps.b(:));
hold(axh,'on')
hold(axh2,'on')
% for nbd = 1:xps.Nbd
%     for nbu = 1:xps.Nbu
%         col = abs(xps.bu_unique(nbu,:));
%         ind = all([nbd==xps.bd_ind nbu==xps.bu_ind],2);
% for nsbdu = 1:xps.Nsbdu
%         nbd = xps.sbdu_unique(nsbdu,1) + 1;
%         seriestemp = xps.sbdu_unique(nsbdu,1);
%         bdtemp = xps.sbdu_unique(nsbdu,2);
%         utemp = xps.sbdu_unique(nsbdu,3:5);

m_v = [m_dti_euler, zeros(1,3*xps.Nbdu)];

for nbdu = 1:xps.Nbdu
        nbd = xps.bdu_unique(nbdu,1) + 1;
        bdtemp = xps.bdu_unique(nbdu,1);
        utemp = xps.bdu_unique(nbdu,2:4);
        col = abs(utemp);
        ind = nbdu==xps.bdu_ind;
%         btemp = [xps.b(xps.b0_logical); xps.b(ind)];
%         Stemp = [S(xps.b0_logical); S(ind)];
%         S_fittemp = [S_fit(xps.b0_logical); S_fit(ind)];
        btemp1 = xps.b(ind);
        Stemp1 = S(ind);
        S_fittemp1 = S_fit(ind);
        seriestemp1 = xps.s_ind(ind);
        [btemp1_unique,~,btemp1_ind] = uniquetol(btemp1,0.01);
        [seriestemp1_unique,~,seriestemp1_ind] = uniquetol(seriestemp1,0.01);
        
        Nrep = numel(btemp1)/numel(btemp1_unique);
        Nseries = numel(seriestemp1_unique);
        
        Nrep = 1; %One fit to all repeated measurements with same b_delta and u
        m_temp = [0 0 0];
        for nrep = 1:Nrep
            if Nrep > 1
                indrep = nrep==seriestemp1_ind;
                btemp = btemp1(indrep);
                Stemp = Stemp1(indrep);
                S_fittemp = S_fittemp1(indrep);
            else
                btemp = btemp1;
                Stemp = Stemp1;
                S_fittemp = S_fittemp1;
            end
            [~,sort_ind] = sort(btemp);
            ph = semilogy(axh2,btemp(sort_ind)+bmax*(nbd-1),Stemp(sort_ind),'o');
            set(ph,'Color',col)
%             ph = semilogy(axh2,btemp(sort_ind)+bmax*(nbd-1),S_fittemp(sort_ind),'-');
%             set(ph,'Color',col)
            %pause(.1)


            if all([sum(ind) > 2 all(btemp>0)])
                weight = ones(size(Stemp));
                t_temp = lsqcurvefit(@my_1d_fit2data, t_guess, btemp, Stemp.*weight,t_lb,t_ub,options);
                s_gamma = my_1d_fit2data(t_temp, btemp);
                m_temp = t2m(t_temp);

                wthresh = opt.color_gamma.weight_sthresh;
                sthresh = opt.color_gamma.weight_wthresh*m_temp(1);
                weight = .5*(erf(wthresh*(s_gamma - sthresh)/sthresh)+1);

                t_temp = lsqcurvefit(@my_1d_fit2data, t_guess, btemp, Stemp.*weight,t_lb,t_ub,options);
                weight = ones(size(Stemp));
                s_gamma = my_1d_fit2data(t_temp, btemp);
                m_temp = t2m(t_temp);

                ph = semilogy(axh2,btemp(sort_ind)+bmax*(nbd-1),s_gamma(sort_ind),'-');
                set(ph,'Color',col)

                %[seriestemp bdtemp utemp]
%                 p2costheta = .5*(3*sum(utemp.*dt.lambda33vec)^2 - 1);
%                 xtemp = 2/3*(p2costheta+.5) + bdtemp;
%                 ph = plot(axh,xtemp,m_temp(2)/MD,'o',xtemp,m_temp(3)/MD^2 + 2,'x');
%                 set(ph,'Color',col)
    %             btemp
    %             numel(btemp)
    %             pause

            end
        end
        m_v(1,(1:3)+3*(nbdu-1)+7) = m_temp;
        
%    end
end

set(axh2,'YScale','log','YLim',[max([.8*min(S(:)); .01*Smax]) 1.1*Smax],'XLim',bmax*[0 2],'XTick',round(bmax*[0 .5 1]/1e8)*1e8)
ylabel(axh2,'signal')
xlabel(axh2,'b / m^-^2s')

m = m_v(1,1:7);
s0          = m(1);
lambdax     = m(2);
lambday     = m(3);
lambdaz     = m(4);
euler_alpha = angle(exp(1i*m(5)));
euler_beta  = angle(exp(1i*m(6)));
euler_gamma = angle(exp(1i*m(7)));

[rotmat,rotmatinv] = tm_euler_angles2rotmat(euler_alpha,euler_beta,euler_gamma);
dt_lambda = diag([lambdax, lambday, lambdaz]);
dt3x3 = rotmat*dt_lambda*rotmatinv;
dt = tm_3x3_to_tpars(dt3x3);

for nbdu = 1:xps.Nbdu
        nbd = xps.bdu_unique(nbdu,1) + 1;
        bdtemp = xps.bdu_unique(nbdu,1);
        utemp = xps.bdu_unique(nbdu,2:4);
        col = abs(utemp);
        m_temp = m_v(1,(1:3)+3*(nbdu-1)+7);
        
        if sum(m_temp)>0
            p2costheta = .5*(3*sum(utemp.*dt.lambda33vec)^2 - 1);
            xtemp = 2/3*(p2costheta+.5) + bdtemp;
            ph = plot(axh,xtemp,m_temp(2)/dt.iso,'o',xtemp,m_temp(3)/dt.iso^2,'x');
            set(ph,'Color',col)
        end
end

plot(axh,[0 1],1*[1 1],'k-')
plot(axh,[1; 2],[dt.lambda11; dt.lambda33]./dt.iso,'k-',[1; 2],[dt.lambda22; dt.lambda33]./dt.iso,'k-')
axis(axh,'tight')
ylim = get(axh,'YLim');
ylim = (ylim(2)-ylim(1))*[-.1 .1]+[ylim(1) ylim(2)];
ylim = ylim(2)*[-.1 1.1];
set(axh,'XLim',[-.1 2.1],'YLim',ylim,'XTick',[0 .5 1])
xlabel(axh,'P_2(cos\theta)')
ylabel(axh,'o D/MD, x \mu_2/MD^2')
title(axh,['FA=' num2str(dt.fa,2) ' u=[' num2str(dt.lambda33vec(1),2) ',' num2str(dt.lambda33vec(2),2) ',' num2str(dt.lambda33vec(3),2) ']'])

set([axh; axh2],'TickDir','out','TickLength',.02*[1 1])

if opt.color_gamma.do_weight

    wthresh = opt.color_gamma.weight_wthresh;
    sthresh = m(1)*opt.color_gamma.weight_sthresh;
    semilogy(axh2, [0; bmax*max(xps.bd_ind)],sthresh*[1; 1],'k--')
%     hold(axh,'on')
%     plot(axh, [0; bmax*max(xps.bd_ind)],sthresh*[1; 1],'k--')

end

end