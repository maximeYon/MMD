clear all

wd = pwd;

DataDir = '/Users/daniel/Dropbox/NMRdata/AVII500/Sarah';
ExpNam = {'Sofia_20180919'}; expno = [5:24];
ExpNam = {'Sofia_20180712'}; expno = [2:37]; expno = 2;
ExpNam = {'YeastForPaper'}; expno = [3 11]; %expno = 11;

for nexp = 1:length(expno)
    data_path = fullfile(DataDir,ExpNam{1},num2str(expno(nexp)));
    
    load(fullfile(data_path,'BSdat'),'BSdat')
    
    
    for peakno = 1:numel(BSdat)
        
%         BSdat{peakno} = struct('m_fit_array',m_fit_array,'chisq',chisq,'s',s,'s_fit',s_fit,'xps',xps);
        m_fit_array = BSdat{peakno}.m_fit_array;
        chisq = BSdat{peakno}.chisq;
        s = BSdat{peakno}.s;
        s_fit = BSdat{peakno}.s_fit;
        xps = BSdat{peakno}.xps;
        
        whiskw = 1.5;
        
        x = m_fit_array(:,9);
        
        med = median(x);
        quant = quantile(x,[.25 .75]);
        q25 = quant(1); q75 = quant(2);
        whisklow = q25 - whiskw*(q75 - q25);
        whiskhigh = q75 + whiskw*(q75 - q25);
        outlier_ind = any([x < whisklow x > whiskhigh],2);
        %figure(1), clf, hist(x(~outlier_ind))
       
        chisq_med = median(chisq);
        chisq_quantile = quantile(chisq,[.25 .75]);
        chisq_q25 = chisq_quantile(1); chisq_q75 = chisq_quantile(2);
        whiskw = 1;
        chisq_whisklow = chisq_q25 - whiskw*(chisq_q75 - chisq_q25);
        chisq_whiskhigh = chisq_q75 + whiskw*(chisq_q75 - chisq_q25);
        %fitOK = all([chisq > chisq_whisklow chisq < chisq_whiskhigh],2);
        fitOK = all([chisq < chisq_whiskhigh ~outlier_ind],2);

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
%%
        s_fit = fexsyr1r2_fit2data(m_fit_array(1,:), xps);

        fh2 = figure(2); clf
        %semilogy(1:xps.n,s,'o',1:xps.n,s_fit,'k-')
        semilogy(1:xps.n,s,'bo')
        xlabel('index'), ylabel('signal')
        set(gca,'Box','off','TickDir','out')
        hold on
        for nbs = 1:numel(fitOK)
            s_fit = fexsyr1r2_fit2data(m_fit_array(nbs,:), xps);
            if fitOK(nbs)
                semilogy(0:xps.n,[sum(m_fit_array(nbs,1:2)); s_fit],'k-')
            else
                %semilogy(0:xps.n,[sum(m_fit_array(nbs,1:2)); s_fit],'r-')
            end
        end
        semilogy(1:xps.n,s,'bo')
       %%
        set(gca,'YLim',sum(m_fit_med(1:2))*[.0002 2], 'XLim',xps.n*[-.1 1.1])
        

        set(fh2, 'PaperUnits','centimeters','PaperPosition', figwidth*[0 0 1 1/aspect],'PaperSize', figwidth*[1 1/aspect]);
        print(fullfile(data_path,['signalfit_peak' num2str(peakno)]),'-dpdf')

    end    
end