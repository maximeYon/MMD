function dtod_nx6 = dtod_nodes2nx6(dtod_nodes,xps)

[n,dpar,dperp,theta,phi,d0,rpar,rperp] = dtod_nodes2par(dtod_nodes);

domega = repmat(xps.domega(1),[1 size(xps.btomega,2)/6]);
omega = cumsum([0 domega(1,1:(end-1))],2);
asize = [1 size(omega,2)];
dtodpars.omega = repmat(omega,[n 1]);
dtodpars.dpar = repmat(dpar,asize);
dtodpars.dperp = repmat(dperp,asize);
dtodpars.theta = repmat(theta,asize);
dtodpars.phi = repmat(phi,asize);
dtodpars.d0 = repmat(d0,asize);
dtodpars.rpar = repmat(rpar,asize);
dtodpars.rperp = repmat(rperp,asize);

dtodpars.dparo = dtodpars.d0 - (dtodpars.d0 - dtodpars.dpar)./(1 + dtodpars.omega.^2./dtodpars.rpar.^2);
dtodpars.dperpo = dtodpars.d0 - (dtodpars.d0 - dtodpars.dperp)./(1 + dtodpars.omega.^2./dtodpars.rperp.^2);

clear dtod
[dtod.xx,dtod.yy,dtod.zz,dtod.xy,dtod.xz,dtod.yz] = dtd_pars2elements(dtodpars.dparo,dtodpars.dperpo,dtodpars.theta,dtodpars.phi);
dtod_nx6 = [dtod.xx dtod.yy dtod.zz sqrt(2)*[dtod.xy dtod.xz dtod.yz]];

% lw = 2;
% lw_v = [lw .75*lw .5*lw .33*lw .33*lw .33*lw lw];
% col_a = [1 0 0
% 0 .7 0
% 0 0 1
% .8 .8 0
% .8 0 .8
% 0 .8 .8
% 0 0 0];
% 
% 
% figure(1), clf
% axh = axes('position',[.1 .1 .8 .8]);
% hold(axh,'on')
% ph_v = [];
% elements = {'xx','yy','zz','xy','xz','yz'};
% for nelement = 1:numel(elements)
%     ph = plot(omega/2/pi,dtod.(elements{nelement})(1,:));
%     ph_v = cat(1,ph_v,ph);
% end
% set(axh,'YLim',max(1*dtodpars.d0(:))*[-.1 1.1],'XLim',max(omega(:)/2/pi)*[-1 1])
% %axis off
% for nline = 1:numel(ph_v)
%     set(ph_v(nline),'LineWidth',lw_v(nline),'Color',col_a(nline,:))
% end

return
xcos = cos(phi).*sin(theta);
ycos = sin(phi).*sin(theta);
zcos = cos(theta);

trace = par + 2*perp;
delta = (par - perp)./trace;

xx = trace/3.*(1 + delta.*(3*xcos.*xcos - 1));
xy = trace/3.*(0 + delta.*(3*xcos.*ycos - 0));
xz = trace/3.*(0 + delta.*(3*xcos.*zcos - 0));
yy = trace/3.*(1 + delta.*(3*ycos.*ycos - 1));
yz = trace/3.*(0 + delta.*(3*ycos.*zcos - 0));
zz = trace/3.*(1 + delta.*(3*zcos.*zcos - 1));

dtod_nx6 = [xx yy zz sqrt(2)*[xy xz yz]];
