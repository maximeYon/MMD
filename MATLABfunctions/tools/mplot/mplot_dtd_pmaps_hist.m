function mplot_dtd_pmaps_hist(voxels, paths, opt)
% function mplot_dtd_pmaps_hist(voxels, paths, opt)

% Font sizes etc
figscale = 2;
figwidth = figscale*17.78;
fs = figscale*6;
lw = figscale*1;
aspect = 1.618;
ms_max = figscale*20;
ms = figscale*3;


param_dist = {'dpar','dperp','theta','phi','w'};
param_map = {'s0','mdiso','msddelta','vdiso','vsddelta','cvdisosddelta'};
param = {param_dist{:}, param_map{:}, 'chisq'};
for nparam = 1:numel(param)
    eval([param{nparam} ' = [];']);
end
for nparam = 1:numel(param_map)
    eval([param_map{nparam} ' = [];']);
end
for nbin = 1:numel(opt.dtd.bin_disomax)
    bin{nbin}.no = nbin;
    for nparam = 1:numel(param_map)
        eval(['bin{nbin}.' param_map{nparam} ' = [];']);
    end
end

Nbins = numel(bin);
Nparams = numel(param_map);

bsno = msf_getdirno(paths.bs_path);
nn = 0;
for nbs = bsno
    mfs_fn   = fullfile(paths.bs_path, num2str(nbs), 'mfs.mat');
    dps_fn   = fullfile(paths.bs_path, num2str(nbs), 'dps.mat');
    chisq_fn   = fullfile(paths.bs_path, num2str(nbs), 'chisq.mat');
    if exist(mfs_fn,'file')==2
        temp = load(mfs_fn); mfs = temp.mfs; 
        m = single(mfs.m);
        %m = mfs.m;
        sz = size(m);
        n = m(:,:,:,1);
        nn_temp = (sz(4)-1)/6;

        ind = false(sz(4),1);
        ind(2:5:end) = 1;

        for nparam = 1:numel(param_dist)
            eval([param_dist{nparam} ' = cat(4,' param_dist{nparam} ',m(:,:,:,circshift(ind,' num2str(nparam-1) ',1)));']);
        end
        nn = nn + nn_temp;
        
        mdm_fit2param(@dtd_4d_fit2param, mfs_fn, dps_fn, opt);
        temp = load(dps_fn); dps = temp.dps; 
        delete(dps_fn);
        for nparam = 1:numel(param_map)
            eval([param_map{nparam} ' = cat(4,' param_map{nparam} ',dps.' param_map{nparam} ');']);
        end
        for nbin = 1:numel(dps.bin)
            for nparam = 1:numel(param_map)
                eval(['bin{nbin}.' param_map{nparam} ' = cat(4,bin{nbin}.' param_map{nparam} ',dps.bin{nbin}.' param_map{nparam} ');']);
            end
        end
    end
    if exist(chisq_fn,'file')==2
        temp = load(chisq_fn); 
        chisq = cat(4,chisq,temp.chisq);
    end
end

for nparam = 1:numel(param_map)
    eval(['dps_bsmean.' param_map{nparam} ' = mean(' param_map{nparam} ', 4);']);
    for nbin = 1:numel(dps.bin)
        eval(['dps_bsmean.bin{nbin}.' param_map{nparam} ' = mean(bin{nbin}.' param_map{nparam} ', 4);']);
    end
end

% Figure 5

figure(5), clf
axh_v = dtd_dpsbins2maps(gcf,dps_bsmean,opt.dtd.plim);
papersize = 3*[Nparams Nbins+1];
set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
fig_fn = fullfile(paths.figs,'Figure5_pmaps');
eval(['print ' fig_fn ' -dpdf -loose'])

% Figure 6

pix_i = voxels.i;
pix_j = voxels.j;
pix_k = voxels.k;

% pixel color
pix_col = [0 0.8 0
         1 .3 .3
         0 0 1];

nbin_tot = numel(bin)+1;
for nparam = 1:numel(param_map)
    eval(['bin{nbin_tot}.' param_map{nparam} ' = ' param_map{nparam} ';']);
end

bin = bin(fliplr(circshift(1:nbin_tot,[1,1])));

s0max = max(dps_bsmean.s0(:));
s0thresh = 0.01;

left = .05; bottom = .07; dleft = (1-left)/Nparams; dheight = (1-bottom)/(Nbins+1);
width = .8*dleft; height = .9*dheight;

Nhistbins = 100;

plim = opt.dtd.plim;
plim.s0 = s0max*[0 1];
pscale.s0 = s0max;
pscale.mdiso = 1e-9;
pscale.msddelta = 1;
pscale.vdiso = 1e-18;
pscale.vsddelta = 1;
pscale.cvdisosddelta = 1e-9;

figure(6), clf

for nparam = 1:Nparams
    for nbin = 1:(Nbins+1)
        dps_temp = bin{nbin};
        eval(['mapdat = dps_temp.' param_map{nparam} ';']);
        eval(['xxlimit = plim.' param_map{nparam} ';'])
        eval(['xxscale = pscale.' param_map{nparam} ';'])
        mapdat = mapdat./xxscale;
        xlim = xxlimit./xxscale;

        axh = axes('position',[left+(nparam-1)*dleft (nbin-1)*dheight+bottom width height]);
        edges = linspace(min(xlim),max(xlim),Nhistbins);
        for npix = 1:3
            ni = pix_i(npix); nj = pix_j(npix); nk = pix_k(npix);
            col = pix_col(npix,:);
            histdat = squeeze(mapdat(ni,nj,nk,:));
            histdats0 = squeeze(dps_temp.s0(ni,nj,nk,:));
            if ~strcmp(param_map{nparam},'s0')
                histdat = histdat(histdats0>s0thresh*s0max);
            end
            h = histogram(histdat,edges,'DisplayStyle','stairs');
            set(h,'EdgeColor',col,'LineWidth',.25*lw*(4-npix))
            set(h,'LineWidth',lw)
            hold on
        end
        set(axh,'XLim',xlim,'YLim',numel(bsno)*[-.1 1.1],'Box','off','TickDir','out','TickLength',.03*[1 1],'LineWidth',lw,'FontSize',fs)
        if nbin > 1
            set(axh,'XTickLabel',[])
        end
        if nparam > 1
            set(axh,'YTickLabel',[])
        end
    end
end

papersize = figwidth*[1 1/aspect];
set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
fig_fn = fullfile(paths.figs,'Figure6_pixelshistograms');
eval(['print ' fig_fn ' -dpdf -loose'])
