function clim = mplot_clim(clim)
% function clim = mplot_clim(clim)
%
% Specifies default color bar limits 

if (nargin < 1), clim.present = 1; end

clim = msf_ensure_field(clim, 's0', .8*[0 1]);
clim = msf_ensure_field(clim, 'mdiso',  3e-9*[0 1]);
clim = msf_ensure_field(clim, 'msddelta', 1*[0 1]);
clim = msf_ensure_field(clim, 'mr1', 10*[0 1]);
clim = msf_ensure_field(clim, 'mr2', 150*[0 1]);

clim = msf_ensure_field(clim, 'vdiso', .5e-18*[0 1]);
clim = msf_ensure_field(clim, 'vsddelta', .15*[0 1]);
clim = msf_ensure_field(clim, 'vr1', 10*[0 1]);
clim = msf_ensure_field(clim, 'vr2', 1000*[0 1]);

clim = msf_ensure_field(clim, 'mask_threshold', eps);

clim = msf_ensure_field(clim, 'cvdisosddelta', sqrt(max(clim.vdiso)*max(clim.vsddelta))*[-1 1]);
clim = msf_ensure_field(clim, 'cvdisor1', sqrt(max(clim.vdiso)*max(clim.vr1))*[-1 1]);
clim = msf_ensure_field(clim, 'cvdisor2', sqrt(max(clim.vdiso)*max(clim.vr2))*[-1 1]);
clim = msf_ensure_field(clim, 'cvsddeltar1', sqrt(max(clim.vsddelta)*max(clim.vr1))*[-1 1]);
clim = msf_ensure_field(clim, 'cvsddeltar2', sqrt(max(clim.vsddelta)*max(clim.vr2))*[-1 1]);
clim = msf_ensure_field(clim, 'cvr1r2', sqrt(max(clim.vr1)*max(clim.vr2))*[-1 1]);

clim = msf_ensure_field(clim, 'dmdisodnu', 1e-3*max(clim.mdiso)*[-1 1]);
clim = msf_ensure_field(clim, 'dmsddeltadnu', 1e-3*max(clim.msddelta)*[-1 1]);
clim = msf_ensure_field(clim, 'dvdisodnu', 1e-3*max(clim.vdiso)*[-1 1]);
clim = msf_ensure_field(clim, 'dvsddeltadnu', 1e-3*max(clim.vsddelta)*[-1 1]);
clim = msf_ensure_field(clim, 'dcvdisosddeltadnu', 1e-3*max(clim.cvdisosddelta)*[-1 1]);

if (~isfield(clim, 'bin')), clim.bin.disomin = [0 0 2]*1e-9; end
clim.bin = msf_ensure_field(clim.bin, 'disomin', [0 0 2]*1e-9);
clim.bin = msf_ensure_field(clim.bin, 'disomax', [2 2 5]*1e-9);
clim.bin = msf_ensure_field(clim.bin, 'sddeltamin', [.25 0 0]);
clim.bin = msf_ensure_field(clim.bin, 'sddeltamax', [1 .25 1]);
clim.bin = msf_ensure_field(clim.bin, 'dratiomin', ones(size(clim.bin.disomin))*eps);
clim.bin = msf_ensure_field(clim.bin, 'dratiomax', ones(size(clim.bin.disomin))/eps);
clim.bin = msf_ensure_field(clim.bin, 'r1min', 0*ones(size(clim.bin.disomin)));
clim.bin = msf_ensure_field(clim.bin, 'r1max', 1/eps*ones(size(clim.bin.disomin)));
clim.bin = msf_ensure_field(clim.bin, 'r2min', 0*ones(size(clim.bin.disomin)));
clim.bin = msf_ensure_field(clim.bin, 'r2max', 1/eps*ones(size(clim.bin.disomin)));



