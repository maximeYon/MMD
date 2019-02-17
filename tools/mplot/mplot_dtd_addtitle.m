function axh = mplot_dtd_addtitle(dps, axh, opt)
% function axh = mplot_dtd_addtitle(dps, axh, opt)
%

switch opt.mplot.dtd_nomenclature
    case 'lasic14'    
        title_str2 = {...
            ['mean "size": MD = ' num2str(dps.mdiso/1e-9, 2) '\cdot10^{-9} m^2/s'];
            ['mean-square "shape": \Delta\mu_2/MD^2 = ' num2str(4/5*dps.msdanison, 2)];
            ['variance "size": \mu_2^{iso}/MD^2 = ' num2str(dps.vdison, 2)]
            [' ',char(181),'FA = ' num2str(dps.ufa, 3)]
        };
    case 'westin16'    
        title_str2 = {...
            ['mean "size": MD = ' num2str(dps.mdiso/1e-9, 2) '\cdot10^{-9} m^2/s'];
            ['mean-square "shape": C_\mu = ' num2str(dps.ufa^2, 2)];
            ['variance "size": C_{MD} = ' num2str(dps.vdison, 2)]
        };
    case 'szczepankiewicz16'    
        title_str2 = {...
            ['mean "size": MD = ' num2str(dps.mdiso/1e-9, 2) '\cdot10^{-9} m^2/s'];
            ['mean-square "shape": MK_a = ' num2str(dps.MKa, 2)];
            ['variance "size": MK_i = ' num2str(dps.MKi, 2)]
        };
    case 'topgaard17'    
        title_str2 = {...
            ['mean "size": E[\itD\rm_{iso}] = ' num2str(dps.mdiso/1e-9, 2) '\cdot10^{-9} m^2/s'];
            ['mean-square "shape": E[(\itD\rm_A-\itD\rm_R)^2]/E[\itD\rm_{iso}]^2 = ' num2str(dps.msdanison, 2)];
            ['variance "size": Var[\itD\rm_{iso}]/E[\itD\rm_{iso}]^2 = ' num2str(dps.vdison, 2)]
        };
    otherwise
        title_str2 = {};
end

title(axh, title_str2,'FontSize',opt.mplot.fs,'FontWeight','normal')

