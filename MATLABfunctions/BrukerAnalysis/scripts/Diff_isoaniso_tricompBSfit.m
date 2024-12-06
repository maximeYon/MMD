clear all

wd = cd;

%DataDir = '/Users/daniel/NMRdata/AVII500/DT';
%DataDir = '/opt/topspin2/data/DT/nmr';
DataDir = '/Users/daniel/Dropbox';

ExpNam = 'Yeast_tripgse'; expno = [73 75 77]; expno = 73; WaterPeak = 1;
signal = 'area'; %signal = 'intensity';
td1start = 2;

MakeBS = 0; %Bootstrap on/off
NBS = 10000; %Number of bootstrap samples
ConfLev = 0.95; %Confidence level

cd(DataDir)
cd(ExpNam)
for nexp = 1:length(expno)
    eval(['load ' num2str(expno(nexp)) '/NMRacqus'])
    eval(['load ' num2str(expno(nexp)) '/FitDat'])
    eval(['load ' num2str(expno(nexp)) '/TricompFitDat'])

    if MakeBS
        if strcmp(signal,'area')        
            S_vector = PP.Apeak(:,WaterPeak);
        else
            S_vector = PP.Ipeak(:,WaterPeak);
        end
        
        Npoints = AcqDat.td1;
        S0 = max(S_vector);

        S0guess = S0*[FitDat.Pout.I0_1 FitDat.Pout.I0_2 FitDat.Pout.I0_3]/FitDat.Pout.I0;
        Dguess = [FitDat.Pout.Diso_1 FitDat.Pout.Diso_2 FitDat.Pout.Diso_3];
        DDeltaguess = [FitDat.Pout.DDelta_1 FitDat.Pout.DDelta_2 FitDat.Pout.DDelta_3]; 
 
        indx_all_array = repmat((1:Npoints)',[1 Npoints]);
        S_array = reshape(S_vector,AcqDat.Nb,AcqDat.Ndir);
        Xin = AcqDat.bmat.b(1:AcqDat.Nb);
        Xin2 = AcqDat.bmat.Delta(1:AcqDat.Nb);        

        Funam = 'fDiffVACSYtricomp';
        Pin = [1*S0guess Dguess DDeltaguess];
        vlb = [.5*S0guess .5*Dguess -.5*ones(1,3)];
        vub = [2*S0guess 2*Dguess 1*ones(1,3)];
        Pnorm = [S0*ones(1,3) Dguess ones(1,3)];

        options = optimset('Display','off');
        Pout_array = zeros(NBS,9);

        tic
        parfor nBS = 1:NBS
            %indx_BS = indx_all;
            indx_BS = td1start - 1 + sort(ceil((Npoints-td1start+1)*rand(Npoints,1)));
            %figure(1), clf, plot(indx_all,indx_BS,'-'), return
            indx_BS_array = repmat(indx_BS',[Npoints 1]);

            PAweight = 1e-3 + sum(indx_all_array == indx_BS_array,2);
            %figure(1), clf, plot(indx_all,PAweight,'-'), return
            PAweight_array = reshape(PAweight,AcqDat.Nb,AcqDat.Ndir);

            S_PA = sum(S_array.*PAweight_array,2)./sum(PAweight_array,2);
            %figure(1), clf, plot(1:AcqDat.Nb*AcqDat.NDeltab,S_PA,'o'), return

            Yin = S_PA; 
            %figure(1), clf, semilogy(Xin,Yin,'o'), return    

            Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); 
            Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,vlb./Pnorm,vub./Pnorm,options,Xin2,Pnorm,Xnorm,Ynorm);
            Yout = feval(Funam,Pout,Xin,Xin2); error = Yin - Yout;
            %figure(1), clf, semilogy(Xin,Yin,'o',Xin,Yout,'x'), return
            Pout_array(nBS,:) = Pout;
        end
        toc

        BSdat.NBS = NBS;
        BSdat.S0_1 = Pout_array(:,1);
        BSdat.S0_2 = Pout_array(:,2);
        BSdat.S0_3 = Pout_array(:,3);
        BSdat.S0 = BSdat.S0_1 + BSdat.S0_2 + BSdat.S0_3;
        BSdat.Diso_1 = Pout_array(:,4);
        BSdat.Diso_2 = Pout_array(:,5);
        BSdat.Diso_3 = Pout_array(:,6);
        BSdat.DDelta_1 = Pout_array(:,7);
        BSdat.DDelta_2 = Pout_array(:,8);
        BSdat.DDelta_3 = Pout_array(:,9);
        BSdat.Dpar_1 = BSdat.Diso_1.*(1 + 2*BSdat.DDelta_1);
        BSdat.Dperp_1 = BSdat.Diso_1.*(1 - BSdat.DDelta_1);
        BSdat.Dpar_2 = BSdat.Diso_2.*(1 + 2*BSdat.DDelta_2);
        BSdat.Dperp_2 = BSdat.Diso_2.*(1 - BSdat.DDelta_2);
        BSdat.Dpar_3 = BSdat.Diso_3.*(1 + 2*BSdat.DDelta_3);
        BSdat.Dperp_3 = BSdat.Diso_3.*(1 - BSdat.DDelta_3);
        
    else
        eval(['load ' num2str(expno(nexp)) '/TricompBSdat'])
    end
%%    
    figure(3), clf
    Nx = 5; Ny = 3;
    ParNam = {'S0','DDelta','Diso','Dpar','Dperp'};
    Ncomp = 3;
    S0tot = 0;
    for ncomp = 1:Ncomp
        eval(['S0tot = S0tot + BSdat.S0_' num2str(ncomp) ';'])
    end
    for ncomp = 1:Ncomp
        for nsub = 1:length(ParNam)
            subplot(Ny,Nx,nsub + length(ParNam)*(ncomp-1))
            if any(strcmp(ParNam{nsub},{'Diso','Dpar','Dperp'}))
                eval(['histdat = log10(BSdat.' ParNam{nsub} '_' num2str(ncomp) ');'])
            elseif any(strcmp(ParNam{nsub},{'S0'}))
                eval(['histdat = BSdat.S0_' num2str(ncomp) './S0tot;'])
            else
                eval(['histdat = BSdat.' ParNam{nsub} '_' num2str(ncomp) ';'])
            end
            histdat = sort(histdat);
            minlim = histdat(floor((1-ConfLev)/2*BSdat.NBS));
            maxlim = histdat(ceil((1-(1-ConfLev)/2)*BSdat.NBS));
            meanval = mean(histdat);
            hist(histdat)       
            %title(['\pm' num2str((maxlim-minlim)/2,2)])
            title([num2str(minlim,4) '-' num2str(maxlim,4)])
            if any(strcmp(ParNam{nsub},{'Diso','Dpar','Dperp'}))
                xlabel(['log(' ParNam{nsub} ')=' num2str(meanval,4)])
            else
                xlabel([ParNam{nsub} '=' num2str(meanval,4)])
            end
            
        end
    end
    
    eval(['save ' num2str(expno(nexp)) '/TricompBSdat BSdat'])
    eval(['print -depsc -loose ' num2str(expno(nexp)) '/TricompBSFig'])
%%        
end

