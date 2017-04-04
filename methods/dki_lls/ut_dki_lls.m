% function fn = ut_dtd_covariance(c_ut)
% % function fn = ut_dtd_covariance(c_ut)
% %
% % Run unit tests on the files in this package
% 
% % if (nargin == 0), fn = 1; return; end
% 
% c_ut = 1;
% 
% switch (c_ut)
%     
%     case 1
%         fn = 'dtd_covariance_1d_fit2data.m';
%         
%         % generate a normally distributed collection of isotropic dffusion 
%         % tensors
%         
%         d_iso = randn(1,1e4) * 0.1e-9 + 1.0e-9;
%         dt = zeros(numel(d_iso), 6);
%         for c = 1:numel(d_iso)
%             dt(c,:) = tm_1x3_to_1x6(1.3 * d_iso(c), d_iso(c), [1 0 0]);
%         end
%         
%         % compute diffusion tensor covariance
%         dtd_cov = mean(tm_1x6_to_1x21(dt),1) - tm_1x6_to_1x21(mean(dt,1));
%         dtd_cov = dtd_cov * 1e18; % convert to natural units (um2/ms)
%         
%         % setup an experiment with LTE and PTE tensors
%         bt = cat(1, ...
%             tm_1x3_to_1x6(eps, 0.0, uvec_tricosa), ...
%             tm_1x3_to_1x6(0.2, 0.0, uvec_tricosa), ...
%             tm_1x3_to_1x6(0.4, 0.0, uvec_tricosa), ...
%             tm_1x3_to_1x6(eps, 0.0, uvec_tricosa), ...
%             tm_1x3_to_1x6(0.0, 0.1, uvec_tricosa), ...
%             tm_1x3_to_1x6(0.0, 0.2, uvec_tricosa)) * 1e9; 
%         
%         xps = mdm_xps_from_bt(bt);
%         
%         xps.u = cat(1, uvec_tricosa, uvec_tricosa, uvec_tricosa, ...
%             uvec_tricosa, uvec_tricosa, uvec_tricosa);
%         
%         % generate the signal
%         s = mean(exp(-xps.bt * dt'), 2);
%         
%         % fit the model
%         opt = dtd_covariance_opt;
%         m = dtd_covariance_1d_data2fit(s, xps, opt);
%         
%         dtd_cov_est = m(8:end)' * 1e18; % convert to um2/ms
%         
%         % check that the estimated cov looks as expected
%         if (numel(dtd_cov_est) ~= numel(dtd_cov))
%             error('%s, ut_dtd_covariance test %i, dtd cov numel wrong', fn, c_ut);
%         end            
%         
%         if (any( abs(dtd_cov_est - dtd_cov) > 1e-4))
%             error('%s, ut_dtd_covariance test %i, dtd cov estimate wrong', fn, c_ut);
%         end
%         
%         
%              
%          
%         
% end
