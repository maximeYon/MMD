function [dtor1r2d_nx6,r1,r2] = dtor1r2d_nodes2nx6r1r2(dtor1r2d_nodes,xps)

[n,dpar,dperp,theta,phi,d0,rpar,rperp,r1,r2] = dtor1r2d_nodes2par(dtor1r2d_nodes);

domega = repmat(xps.domega(1),[1 size(xps.btomega,2)/6]);
omega = cumsum([0 domega(1,1:(end-1))],2);
asize = [1 size(omega,2)];
dtor1r2dpars.omega = repmat(omega,[n 1]);
dtor1r2dpars.dpar = repmat(dpar,asize);
dtor1r2dpars.dperp = repmat(dperp,asize);
dtor1r2dpars.theta = repmat(theta,asize);
dtor1r2dpars.phi = repmat(phi,asize);
dtor1r2dpars.d0 = repmat(d0,asize);
dtor1r2dpars.rpar = repmat(rpar,asize);
dtor1r2dpars.rperp = repmat(rperp,asize);

dtor1r2dpars.dparo = dtor1r2dpars.d0 - (dtor1r2dpars.d0 - dtor1r2dpars.dpar)./(1 + dtor1r2dpars.omega.^2./dtor1r2dpars.rpar.^2);
dtor1r2dpars.dperpo = dtor1r2dpars.d0 - (dtor1r2dpars.d0 - dtor1r2dpars.dperp)./(1 + dtor1r2dpars.omega.^2./dtor1r2dpars.rperp.^2);

clear dtor1r2d
[dtor1r2d.xx,dtor1r2d.yy,dtor1r2d.zz,dtor1r2d.xy,dtor1r2d.xz,dtor1r2d.yz] = dtd_pars2elements(dtor1r2dpars.dparo,dtor1r2dpars.dperpo,dtor1r2dpars.theta,dtor1r2dpars.phi);
dtor1r2d_nx6 = [dtor1r2d.xx dtor1r2d.yy dtor1r2d.zz sqrt(2)*[dtor1r2d.xy dtor1r2d.xz dtor1r2d.yz]];

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
%     ph = plot(omega/2/pi,dtor1r2d.(elements{nelement})(1,:));
%     ph_v = cat(1,ph_v,ph);
% end
% set(axh,'YLim',max(1*dtor1r2dpars.d0(:))*[-.1 1.1],'XLim',max(omega(:)/2/pi)*[-1 1])
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

dtor1r2d_nx6 = [xx yy zz sqrt(2)*[xy xz yz]];
