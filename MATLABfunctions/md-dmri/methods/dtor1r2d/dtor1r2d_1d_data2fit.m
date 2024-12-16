function m = dtor1r2d_1d_data2fit(signal, xps, opt, ind)
% function m = dtor1r2d_1d_data2fit(signal, xps, opt, ind)
%
% Size-shape-orientation-frequency diffusion tensor distribution

if (nargin < 4), ind = ones(size(signal)) > 0; end

stemp = signal(ind);
stemp = abs(stemp);
xpstemp.n = numel(stemp);
xpstemp.btomega = double(xps.btomega(ind,:));
xpstemp.domega = xps.domega(ind,:);
xpstemp.te = xps.te(ind,:);
xpstemp.tr = xps.tr(ind,:);

dtor1r2d = dtor1r2d_proliferation(stemp, xpstemp, opt);
%dtor1r2d
%pause
% dpar = [.01e-9; 1e-9];
% dperp = [.01e-9; 1e-9];
% theta = [pi/5; 0];
% phi = [pi/8; 0];
% d0 = [1e-9; 1e-9];
% rpar =  [500; 1];
% rperp = [500; 1];
% w = 7e9*[1; 1];
% dtor1r2d = dtor1r2d_par2dist(dpar,dperp,theta,phi,d0,rpar,rperp,w);

dtor1r2d = dtor1r2d_darwin(stemp, xpstemp, dtor1r2d, opt);

m = dtor1r2d_dtor1r2d2m(dtor1r2d,opt);
%size(m)

if (opt.dtor1r2d.do_plot)
    figure(1), clf
    signal_fit = dtor1r2d_1d_fit2data(m, xpstemp);
    %[~,s_ind] = sort(signal_fit,'descend');
    %semilogy(xps.b,signal,'.',xps.b,signal_fit,'o',xps.b,m(1)*weight,'x');
    plot(1:xpstemp.n,stemp,'o',1:xpstemp.n,signal_fit,'x');
%     plot(ind,signal(ind),'o',ind,signal_fit(ind),'x');
%    plot(te,stemp,'o',te,signal_fit(ind),'x');
%    plot(sum(bt_mx6(:,1:3),2),stemp,'o',sum(bt_mx6(:,1:3),2),signal_fit(ind),'x');
    %set(gca,'YLim',m(1)*[.01 1.2])
    pause(0.05);
end