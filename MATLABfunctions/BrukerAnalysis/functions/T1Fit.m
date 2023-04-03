clear FitDat
FitDat.fitpoints = fitpoints;
FitDat.signal = signal;
for npeak = 1:Npeaks
    Xin = Xdat(fitpoints);
    Yin = PP.Ipeak(fitpoints,npeak);
    if strcmp(signal,'area')==1
        Yin = PP.Apeak(fitpoints,npeak);
    end
    %figure(1), clf, semilogy(Xin,Yin,'o'), return

    Pin = [max(Yin)   1/mean(Xin) 1.5]; Funam = 'fIR';

    Pout = Pin; Ynorm = mean(Yin); Xnorm = mean(Xin); Pnorm = Pin; 
    Pout = Pnorm.*lsqcurvefit(Funam,Pin./Pnorm,Xin/Xnorm,Yin/Ynorm,zeros(size(Pin)),[],[],Pnorm,Xnorm,Ynorm);
    Yout = feval(Funam,Pout,Xin,ones(size(Pin)),1,1); error = Yin - Yout;

    if PlotInterm
        figure(1), clf
        axes('position',[.1 .27 .8 .65])
        semilogy(Xin,Yin,'o',Xin,Yout)
        title(['expno=' num2str(expno(nexp)) '  peak=' num2str(npeak) '  '  Funam  '   Pout= ' num2str(Pout,3) ])
        axis([min(Xin) max(Xin) max(Yout)*[1e-3 1.2] ])
        xlabel('time'), ylabel('intensity')
        axes('position',[.1 .05 .8 .1])
        plot(Xin,error,'o'), grid
        ylabel('residual')
        pause(1)
        %return
    end

    FitDat.Xin(:,npeak) = Xin;
    FitDat.Yin(:,npeak) = Yin;
    FitDat.Yout(:,npeak) = Yout;
    FitDat.Y0(:,npeak) = Pout(1);
    FitDat.R(:,npeak) = Pout(2);
end 
