function mfs_fn = dti_nls_4d_data2fit(s, o, opt)
% function mfs_fn = dti_nls_4d_data2fit(s, o, opt)
%
% Loops over a 4D volume to produce fit parameters with the DTI nls model
%
% Input: 
%
% s   - data structure
%          s.nii_fn - full path to nifti filename with data
%          s.xps    - experimental parameter structure
%
% o   - output folder   
% 
% opt - options structure, optional
%
% Output:
%
% mfs_fn - path to .mat file with the ModelFitStructure (mfs)

if (nargin < 3), opt = []; end

opt = dti_nls_opt(opt);

    function m = my_core(signal, do_fit)
        
        if (do_fit)
            m = dti_nls_1d_data2fit(signal, s.xps, opt);
        else
            m = zeros(1,7);
        end
        
        if (opt.dti_nls.do_plot)
            signal_fit = dti_nls_1d_fit2data(m, s.xps);
            x = (1:numel(signal))';
            plot(x,signal,'.',x,signal_fit,'o');
            pause(0.05);
        end
        
    end

% Verify the xps
dti_nls_check_xps(s.xps);

% Loop over the volume and fit the model
mfs_fn = mio_fit_model(@my_core, s, o, opt);

end