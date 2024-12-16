function mplot_globalstats_perbin(method, dps, fig_fn, clim, opt)
% mplot_globalstats_perbin(method, dps, fig_fn, clim, opt)

parnams_struct.means = {};
parnams_struct.vars = {};
parnams_struct.covs = {};
parnams_struct.dnus = {};
parnams_struct.dnus = {};
if strcmp(method,'dtor1r2d')
    parnams_struct.means = {'mdiso';'msddelta';'mr1';'mr2';'f'};
    parnams_struct.vars = {'vdiso';'vsddelta';'vr1';'vr2'};
    parnams_struct.covs = {'cvdisosddelta';'cvdisor1';'cvdisor2';'cvsddeltar1';'cvsddeltar2';'cvr1r2'};
    parnams_struct.dnus = {'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu'};
elseif strcmp(method,'dtr1r2d')
    parnams_struct.means = {'mdiso';'msddelta';'mr1';'mr2';'f'};
    parnams_struct.vars = {'vdiso';'vsddelta';'vr1';'vr2'};
    parnams_struct.covs = {'cvdisosddelta';'cvdisor1';'cvdisor2';'cvsddeltar1';'cvsddeltar2';'cvr1r2'};
%     parnams_struct.covs = {'cvdisosddelta';'cvdisor1';'cvdisor2';'cvr1r2'};
elseif strcmp(method,'dtr1d')
    parnams_struct.means = {'mdiso';'msddelta';'mr1';'f'};
    parnams_struct.vars = {'vdiso';'vsddelta';'vr1'};
    parnams_struct.covs = {'cvdisosddelta';'cvdisor1';'cvsddeltar1'};
elseif strcmp(method,'dtr2d')
    parnams_struct.means = {'mdiso';'msddelta';'mr2';'f'};
    parnams_struct.vars = {'vdiso';'vsddelta';'vr2'};
    parnams_struct.covs = {'cvdisosddelta';'cvdisor2';'cvsddeltar2'};
elseif strcmp(method,'dtod')
    parnams_struct.means = {'mdiso';'msddelta';'f'};
    parnams_struct.vars = {'vdiso';'vsddelta'};
    parnams_struct.covs = {'cvdisosddelta'};
    parnams_struct.dnus = {'dmdisodnu';'dmsddeltadnu';'dvdisodnu';'dvsddeltadnu';'dcvdisosddeltadnu'};
elseif strcmp(method,'dtd')
    parnams_struct.means = {'mdiso';'msddelta';'f'};
    parnams_struct.vars = {'vdiso';'vsddelta'};
    parnams_struct.covs = {'cvdisosddelta'};
end

%Plot histograms
figure(3), clf
opt = mplot_opt();
opt.mplot.fs = 7;

statnams = fieldnames(parnams_struct);
Nparam = max([numel(parnams_struct.means) numel(parnams_struct.vars) numel(parnams_struct.covs) numel(parnams_struct.dnus)]);
Nrow = numel(statnams);
%Nrow = 1; % Skip vars and covs
papersize = [17.56 10/3*Nrow];
left0 = .1/papersize(1);
bottom0 = .8/papersize(2);
dleft = (1-left0)/(Nparam+0);
dbottom = 1/Nrow;
width = .9*dleft;
height = (1-Nrow*bottom0)/Nrow;

if ~isfield(dps,'f')
    dps.f = ones(size(dps.s0));
    clim.f = [0 1];
end

for nfield = 1:Nrow
    statnam = statnams{nfield};
    parnams = parnams_struct.(statnam);
    bottom = bottom0+dbottom*(Nrow-nfield);
    for nparam = 1:numel(parnams)    
        parnam = parnams{nparam};
        left = left0+dleft*(nparam-1);
        axh = axes('position',[left bottom width height]);
        [axh,hists] = plot_hist_perbin(axh,parnam);
        xlabel(parnam)
        hists_struct.(parnam) = hists;
    end
end

set(gcf, 'PaperUnits','centimeters','PaperPosition', [0 0 papersize],'PaperSize', papersize); 
set(gcf, 'InvertHardcopy','off','Color',[1 1 1]); 

if ~isempty(fig_fn)
    msf_mkdir(fileparts(fig_fn));
    print(fig_fn,'-loose','-dpdf')
    save(fullfile(fileparts(fig_fn),[method '_hists']),'hists_struct')
end


    function [axh,hists] = plot_hist_perbin(axh,parnam)
        xlim = clim.(parnam) + .1*abs(diff(clim.(parnam)))*[-1 1];
        dist_s.x = linspace(min(xlim),max(xlim),200)';
        dist_s.xsigma = 1*(dist_s.x(2) - dist_s.x(1)); % 3 instead of 1 

        dist_d.w = dps.s0(:);
        dist_d.x = dps.(parnam)(:);

        dist_d.n = numel(dist_d.w);
        dist_s = dist_1d_discrete2smooth(dist_d,dist_s);

        wmax = max(dist_s.w);

        plot(axh,dist_s.x, dist_s.w,'-','LineWidth',2*opt.mplot.lw,'Color',.8*[1 1 1])
        hold(axh,'on')
        
        hists.glob = dist_s;

        color_a = [[1 0 0]; [0 1 0]; [0 0 1]];
        wmax_perbin = zeros(3,1);
        for nbin = 1:3
            dist_d.w = dps.bin{nbin}.s0(:);
            dist_d.x = dps.bin{nbin}.(parnam)(:);
            dist_s = dist_1d_discrete2smooth(dist_d,dist_s);
            plot(axh,dist_s.x, dist_s.w,'-','LineWidth',1.5*(6-nbin)/6*opt.mplot.lw,'Color',color_a(:,nbin))
            
            hists.bin{nbin} = dist_s;
            
            wmax_perbin(nbin,1) = max(dist_s.w);
        end

        axis(axh,'tight')
        set(axh,'Box','off','TickDir','out','TickLength',.02*[1 1],...
            'YTick',[],'YColor','none','XLim',xlim,'YLim',max(wmax_perbin)*[-.1 1.5],'FontSize',opt.mplot.fs,'LineWidth',opt.mplot.lw)
    end

end
