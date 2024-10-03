I = real(I);
Itd1 = real(Itd1);

if exist('td1start') == 0
    Proc.td1start = 1;
else
    Proc.td1start = td1start;  
end
if any(strcmp(NMRacqus.pulprog,{'DT_T1ir','DT_seT1ir'})) == 1
    Proc.td1start = NMRacqu2s.td;
end

if exist('thresh') == 0
    Proc.thresh = .1;
else
    Proc.thresh = thresh;  
end

ntd1 = Proc.td1start;

if exist('si') == 0
    Proc.si = 2*td/2;
else
    Proc.si = si;    
end

Itemp = zeros(Proc.si,1);
indx = round(.01*Proc.si):round(.99*Proc.si);
Itemp(indx,1) = Itd1(indx,ntd1);
Imax = max(Itemp);


if any([FindPeaks exist('ul') == 0]) == 1
    pos = find(Itemp>(Imax*Proc.thresh));
    posn = find(diff(pos)~=1);
    ul=[pos(posn); max(pos)];
    ll=[min(pos); pos(posn+1)];
    Npeaks = length(ul);

    del = [];
    for npeak = 1:Npeaks
        if (ul(npeak)-ll(npeak))<5
            del = [del; npeak];
        end
    end
    ll(del) = [];
    ul(del) = [];
    Npeaks = length(ul);

    del = [];
    for npeak = 1:(Npeaks-1)
        if (ll(npeak+1)-ul(npeak))<5
            del = [del; npeak];
        end
    end
    ll(del+1) = [];
    ul(del) = [];

end
Npeaks = length(ul);

intlim = [ll'; ul'];

if CheckPeaks
    PlotInterm = 1;
end

if PlotInterm
    clear Ipeak peakpos
    for npeak = 1:Npeaks
        peak = Itemp(ll(npeak):ul(npeak),1);
        Ipeak(ntd1,npeak) = max(peak);
        peakpos(ntd1,npeak) = max(find(Itemp==Ipeak(ntd1,npeak)));
    end
    
    figure(3), clf
    if [Npeaks>1 & max(Ipeak(ntd1,:))>4*min(Ipeak(ntd1,:))]
        subplot(2,2,1)
        plot(1:Proc.si,Itemp,peakpos(ntd1,:),Ipeak(ntd1,:),'s',intlim,zeros(size(intlim)),'r-',baslpoints,zeros(size(baslpoints)),'k.',[1 Proc.si],Imax*Proc.thresh*[1 1],'k--')
        ax = [0 Proc.si Imax*[-.1 1.1]];
        xlabel('channel')
        title(['Peak pick   expno=' num2str(expno(nexp)) '  td1start=' num2str(td1start)])
        axis(ax)
        set(gca,'XDir','reverse')
        subplot(2,2,2)
        plot(1:Proc.si,Itemp,peakpos(ntd1,:),Ipeak(ntd1,:),'s',intlim,zeros(size(intlim)),'r-',baslpoints,zeros(size(baslpoints)),'k.',[1 Proc.si],Imax*Proc.thresh*[1 1],'k--')
        xlabel('channel')
        axis([.8*min(ll) 1.1*max(ul) abs(min(Ipeak(ntd1,:)))*[-.1 2]])
        set(gca,'XDir','reverse')
    else
        subplot(2,1,1)
        plot(1:Proc.si,Itemp,peakpos(ntd1,:),Ipeak(ntd1,:),'s',intlim,zeros(size(intlim)),'r-',baslpoints,zeros(size(baslpoints)),'k.',[1 Proc.si],Imax*Proc.thresh*[1 1],'k--')
        ax = [0 Proc.si Imax*[-.1 1.1]];
        xlabel('channel')
        title(['Peak pick   expno=' num2str(expno(nexp)) '  td1start=' num2str(Proc.td1start)])
        axis(ax)
        set(gca,'XDir','reverse')
    end

    pause(0.1)
    if CheckPeaks
        return
    end
end

Ipeak = zeros(td1,Npeaks);
Apeak = zeros(td1,Npeaks);
peakpos = zeros(td1,Npeaks);
peakppm = zeros(td1,Npeaks);
for ntd1 = 1:td1
    Itemp(:,1) = Itd1(:,ntd1);
    

    for npeak = 1:Npeaks
        peak = Itemp(ll(npeak):ul(npeak));
        Ipeak(ntd1,npeak) = max(peak);
        Apeak(ntd1,npeak) = sum(peak);
        peakpos(ntd1,npeak) = max(find(Itemp==Ipeak(ntd1,npeak)));
        peakppm(ntd1,npeak) = ppm(peakpos(ntd1,npeak));
    end
    if PlotInterm
        figure(1), clf
        plot(1:Proc.si,Itemp,peakpos(ntd1,:),Ipeak(ntd1,:),'s',intlim,zeros(size(intlim)),'r-',baslpoints,zeros(size(baslpoints)),'k.',[1 Proc.si],Imax*Proc.thresh*[1 1],'k')
        ax = [0 Proc.si Imax*[-.1 1.1]];
        xlabel('channel')
        title(['Peak pick   expno=' num2str(expno(nexp)) '  td1=' num2str(ntd1)])
        set(gca,'XDir','reverse')
        pause(.1)
    end
end

if PlotInterm
    figure(3)
    subplot(2,2,3), semilogy(Ipeak(Proc.td1start:td1,:),'o')
    title('intensity'), xlabel('td1'), axis tight
    subplot(2,2,4), semilogy(Apeak(Proc.td1start:td1,:),'o')
    title('area'), xlabel('td1'), axis tight
    pause(.1)
end

PP = struct('Ispec',Itd1(:,Proc.td1start),'ppm',ppm,'intlim',intlim,'Npeaks',Npeaks,...
'Ipeak',Ipeak,'Apeak',Apeak,'peakindx',peakpos,'peakppm',peakppm,'td1start',Proc.td1start);

clear Proc