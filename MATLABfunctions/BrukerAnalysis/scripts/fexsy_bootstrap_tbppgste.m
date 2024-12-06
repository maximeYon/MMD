clear all

wd = pwd;

DataDir = '/Users/daniel/Dropbox/NMRdata/AVII500/Sarah';
ExpNam = {'Sofia_20180919'}; expno = [5:24];
ExpNam = {'Sofia_20180712'}; expno = [2:37]; expno = 2;
ExpNam = {'YeastForPaper'}; expno = [3 11];

Nbs = 400;

for nexp = 1:length(expno)
    data_path = fullfile(DataDir,ExpNam{1},num2str(expno(nexp)));

    load(fullfile(data_path,'NMRacqu2s'))
    td1 = NMRacqu2s.td;
    load(fullfile(data_path,'NMRacqus'))

    % Load gyromagnetic ratio gamma
    gamma = mdm_bruker_gamma(NMRacqus);

    % Load max gradient Gmax
    Gmax = mdm_bruker_maxgradient(NMRacqus);

    gnams = {'ax','bx','cx','ay','by','cy','az','bz','cz'};
    for ngnam = 1:numel(gnams)
        gnam = gnams{ngnam};
        fn = fullfile(data_path,['g' gnam]);
        g = mdm_bruker_grad_read(fn);
        ramp.(gnam) = g;
    end

    fid = fopen([data_path '/vdT1']);
    vdT1 = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/vdT2']);
    vdT2 = fscanf(fid,'%f');
    fclose(fid);
    fid = fopen([data_path '/vdTM']);
    vdTM = fscanf(fid,'%f');
    fclose(fid);

    epsilon = NMRacqus.d2;
    tau = 2*NMRacqus.d4 + 1e-6*NMRacqus.p2;
    delta = 2*(NMRacqus.d2 + NMRacqus.d3);
    Delta = 2*NMRacqus.d2 + 2*NMRacqus.d3 + 2*NMRacqus.d4 + 2*NMRacqus.d6 + NMRacqus.d5;            
    DeltaM = 2*NMRacqus.d2 + 2*NMRacqus.d3 + 2*NMRacqus.d4 + 2*NMRacqus.d6 + vdTM;            
    if NMRacqus.l11
        Delta = Delta + NMRacqus.d46 + 2*NMRacqus.d42 + NMRacqus.d43;
        DeltaM = DeltaM + NMRacqus.d46 + 2*NMRacqus.d42 + NMRacqus.d43;
    end
    tdiff = Delta - delta/3 - tau/2 - epsilon/2 - epsilon^2/6/delta + epsilon^3/15/delta^2;
    tdiffM = DeltaM - delta/3 - tau/2 - epsilon/2 - epsilon^2/6/delta + epsilon^3/15/delta^2;

    G1 = Gmax*sqrt((NMRacqus.cnst1/100*ramp.ax).^2+(NMRacqus.cnst2/100*ramp.ay).^2+(NMRacqus.cnst3/100*ramp.az).^2);
    b1 = gamma^2*G1.^2*delta^2*tdiff;
    Gm = Gmax*sqrt((NMRacqus.cnst1/100*ramp.bx).^2+(NMRacqus.cnst2/100*ramp.by).^2+(NMRacqus.cnst3/100*ramp.bz).^2);
    bm = gamma^2*Gm.^2*delta^2.*tdiffM;
    G2 = Gmax*sqrt((NMRacqus.cnst1/100*ramp.cx).^2+(NMRacqus.cnst2/100*ramp.cy).^2+(NMRacqus.cnst3/100*ramp.cz).^2);
    b2 = gamma^2*G2.^2*delta^2*tdiff;

    tT2 = 2*(NMRacqus.d6 + 2*NMRacqus.d2 + NMRacqus.d3 + NMRacqus.d4);
    tT1 = NMRacqus.d5;
    if NMRacqus.l11
        tT1 = tT1 + (NMRacqus.d46 + 2*NMRacqus.d42 + NMRacqus.d43);
    end

    xps.n = td1;
    xps.q1 = sqrt(b1/tdiff);
    xps.q2 = sqrt(b2/tdiff);

    xps.vtau_KR1 = vdTM;
    xps.vtau_KR2 = vdT2 + 2*(NMRacqus.d56 + 2*NMRacqus.d52 + NMRacqus.d53 + NMRacqus.d22);
    xps.tau_KR2 = (tT2-(Delta-tT1)/2)*ones(xps.n,1);
    xps.tau_KR2D = (Delta-tT1)/2*ones(xps.n,1);
    xps.tau_KR1D = tT1*ones(xps.n,1);
    xps.tr = vdT1;

%     figure(1), clf
%     subplot(2,1,1)
%     plot(1:xps.n,b1,'o',1:xps.n,b2,'o',1:xps.n,bm,'o')
%     subplot(2,1,2)
%     semilogy(1:xps.n,xps.vtau_KR1,'o',1:xps.n,xps.vtau_KR2,'o',1:xps.n,xps.tau_KR2,'o',1:xps.n,xps.tau_KR1D,'o',1:xps.n,xps.tr,'o')

    load(fullfile(data_path,'PeakDat'))

    full_xps = xps;
    for peakno = 1:PP.Npeaks
    %for peakno = 2

        opt.fexsyr1r2.present = 1;
        opt.fexsyr1r2 = msf_ensure_field(opt.fexsyr1r2, 'lsq_opts', ...
            optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e3));

        s = PP.Apeak(:,peakno);

        if abs((PP.peakppm(1,peakno)-2.63))<.2 % glycerol peak with J-couplings
            sub_ind = xps.vtau_KR2 < 10e-3;
            s = s(sub_ind);
            xps = mdm_xps_subsample(full_xps,sub_ind); xps.n = numel(xps.q1);
        else
            xps = full_xps;
        end
        
        ind = 1:xps.n;

        ind_start = 1;
        chisq = zeros(Nbs,1);
        m_fit_array = zeros(Nbs,9);
        parfor nbs = 1:Nbs
            ind = (ind_start-0) + round(rand([xps.n-(ind_start-1),1])*(xps.n-(ind_start-0)));
            m_fit = fexsyr1r2_data2fit(s, xps, opt, ind);
            s_fit = fexsyr1r2_fit2data(m_fit, xps);
            chisq(nbs) = sum((s-s_fit).^2)/xps.n;
            m_fit_array(nbs,:) = m_fit;

        %     figure(2), clf
        %     semilogy(1:xps.n,s,'o',1:xps.n,s_fit,'x')
        %     set(gca,'YLim',sum(m_fit(1:2))*[1e-3 1])
        %     pause(.1)
        end


        chisq_med = median(chisq);
        chisq_quantile = quantile(chisq,[.25 .75]);
        chisq_q25 = chisq_quantile(1); chisq_q75 = chisq_quantile(2);
        whiskw = 1;
        chisq_whisklow = chisq_q25 - whiskw*(chisq_q75 - chisq_q25);
        chisq_whiskhigh = chisq_q75 + whiskw*(chisq_q75 - chisq_q25);
        fitOK = all([chisq > chisq_whisklow chisq < chisq_whiskhigh],2);

        % figure(4), clf
        % boxplot(chisq)

        m_fit_array(:,10) = m_fit_array(:,1) + m_fit_array(:,2);
        m_fit_array(:,11) = m_fit_array(:,1)./m_fit_array(:,10);

        xlabel_str = {'S_i'; 'S_e'; 'R_1_i / s^-^1'; 'R_1_e / s^-^1'; 
            'R_2_i / s^-^1'; 'R_2_e / s^-^1'; 'D_i / m^2s^-^1'; 'D_e / m^2s^-^1'; 'k_i_e / s^-^1'; 'S_0'; 'P_i'};

        fh3 = figure(3); clf
        for n = 1:size(m_fit_array,2)
            subplot(3,4,n)
            histogram(m_fit_array(fitOK,n))
            set(gca,'Box','off','TickDir','out')
            xlabel(xlabel_str{n})
        end
        subplot(3,4,12)
        histogram(sqrt(chisq(fitOK)))
        set(gca,'Box','off','TickDir','out')
        xlabel('rms\chi')


        figwidth = 17.78;
        aspect = 1.618;
        set(fh3, 'PaperUnits','centimeters','PaperPosition', figwidth*[0 0 1 1/aspect],'PaperSize', figwidth*[1 1/aspect]);
        print(fullfile(data_path,['histograms_peak' num2str(peakno)]),'-dpdf')

        m_fit_med = median(m_fit_array(fitOK,:),1);
        m_fit_iqr = iqr(m_fit_array,1);
        s_fit = fexsyr1r2_fit2data(m_fit_med, xps);

        fh2 = figure(2); clf
        semilogy(1:xps.n,s,'o',1:xps.n,s_fit,'-')
        set(gca,'YLim',[.1*min(sqrt(chisq)) sum(m_fit_med(1:2))])
        xlabel('index'), ylabel('signal')
        set(gca,'Box','off','TickDir','out')

        set(fh2, 'PaperUnits','centimeters','PaperPosition', figwidth*[0 0 1 1/aspect],'PaperSize', figwidth*[1 1/aspect]);
        print(fullfile(data_path,['signalfit_peak' num2str(peakno)]),'-dpdf')

        BSdat{peakno} = struct('m_fit_array',m_fit_array,'chisq',chisq,'s',s,'s_fit',s_fit,'xps',xps);
    end
    %%
    save(fullfile(data_path,'BSdat'),'BSdat')
end