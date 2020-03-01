function opt = dtd_gamma_opt(opt)
% function opt = dtd_gamma_opt(opt)
%
% Makes sure that all needed fields in the options structure are present

opt.dtd_gamma.present = 1;

x = 'dtd_gamma'; % for shortening 

opt.(x) = msf_ensure_field(opt.(x), 'tmp', 1); 
opt.(x) = msf_ensure_field(opt.(x), 'pa_method', 1); 
opt.(x) = msf_ensure_field(opt.(x), 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off','MaxFunEvals',1e3));
opt.(x) = msf_ensure_field(opt.(x), 'do_plot', 0);
opt.(x) = msf_ensure_field(opt.(x), 'do_weight', 0);
opt.(x) = msf_ensure_field(opt.(x), 'do_pa_weight', 1);
opt.(x) = msf_ensure_field(opt.(x), 'weight_sthresh', .07);
opt.(x) = msf_ensure_field(opt.(x), 'weight_wthresh', 2);
opt.(x) = msf_ensure_field(opt.(x), 'weight_mdthresh', 1e-9);
opt.(x) = msf_ensure_field(opt.(x), 'do_clamping', 0);
opt.(x) = msf_ensure_field(opt.(x), 'fig_maps', ...
    {'s0','MD', 'mu2ison', 'mu2anison', 'uFA'});
opt.(x) = msf_ensure_field(opt.(x), 'fig_prefix', 'dtd_gamma');

% Number of identical fitting repetitions (keeping smallest residual)
% Very expensive
opt.(x) = msf_ensure_field(opt.(x), 'fit_iters', 1);

% Number of random guesses to start from (keeping guess with smallest residual)
% Not expensive
opt.(x) = msf_ensure_field(opt.(x), 'guess_iters', 50);

% Decide if multiple series should be assumed to have different baseline
% signal (s_ind has multiple unique values).
opt.(x) = msf_ensure_field(opt.(x), 'do_multiple_s0', 1);

% Bounds for initial guess and fitting (not including relative signal for
% multiple series. This is done in the data2fit function.
% Note that the first value is S(b=0)/max(signal), i.e., a relative signal.
opt.(x) = msf_ensure_field(opt.(x), 'fit_lb', [  0 1e-12 [1 1]*1e-24 ]);
opt.(x) = msf_ensure_field(opt.(x), 'fit_ub', [ 10 4e-9  [1 1]*5e-18 ]);

% used in dtd_gamma_pipe.m
opt.(x) = msf_ensure_field(opt.(x), 'do_pa', 1);

% used in mdm_fit.m (work needed here to streamline things)
opt.(x) = msf_ensure_field(opt.(x), 'pa_method', 1);
