function mfs_fn = dti_euler_4d_data2fit(s, o, opt)
% function mfs_fn = dti_euler_4d_data2fit(s, o, opt)

if (nargin < 3), opt = []; end

% opt = vasco16_opt(opt);
% opt.vasco16.do_plot = 0;


ind_signal = opt.td1start:s.xps.n;

    function m = my_core(signal, do_fit)
        if (do_fit)
            m = dti_euler_1d_data2fit(signal, s.xps, opt, ind_signal);
        else
            m = zeros(1,7);
        end
        
%         if (opt.vasco16.do_plot) && (do_fit)
%             signal_fit = vasco16_1d_fit2data(m, s.xps);
%             
%             for c = 1:2
%                 subplot(1,2,c);
%                 semilogy(...
%                     s.xps.b(ind{c}) * 1e-6, signal(ind{c}),'x',...
%                     s.xps.b(ind{c}) * 1e-6, signal_fit(ind{c}),'o-'); 
%                 ylim( [0.7 1] * max(signal));
%                 title(t{c});
%             end
%             pause(0.05);
%         end   
    end

% Verify the xps
%vasco16_check_xps(s.xps);

% Loop over the volume and fit the model
mfs_fn = mio_fit_model(@my_core, s, o, opt);

end