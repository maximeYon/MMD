for ntd1 = 1:td1
    Itemp(:,1) = Itd1(:,ntd1);

    Pout = fminsearch('ACME_AutoPhase',[10 -10],[],Itemp,baslpoints,pivotpoint);
    phc0=Pout(1);
    phc1=Pout(2);
    phcorrfun = exp(1i*(phc0*pi/180.*ones(1,si)'+phc1*pi/180*((1:si)'-pivotpoint)./si));
    Itemp = Itemp.*phcorrfun;
    Itemp = Itemp - mean(Itemp(baslpoints));

    Itd1(:,ntd1) = Itemp(:,1);

    if PlotInterm
        figure(1), clf
        plot(1:si,real(Itemp),'k')
        ax = [0 si max(real(Itemp))*[-.1 1.1]];
        xlabel('channel')
        title(['Peak pick   expno=' num2str(expno(nexp)) '  td1=' num2str(ntd1)])
        set(gca,'XDir','reverse')
        pause(.1)
    end
end
