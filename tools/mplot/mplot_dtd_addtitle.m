function axh = mplot_dtd_addtitle(dps, axh, opt)
% function axh = mplot_dtd_addtitle(dps, axh, opt)
%
    switch opt.mplot.terminology
        case 'Lasic14'    
            title_str2 = {...
                ['mean "size": MD = ' num2str(dps.mdiso/1e-9, 2) '\cdot10^{-9} m^2/s'];
                ['mean-square "shape": \Delta\mu_2/MD^2 = ' num2str(4/5*dps.msdanison, 2)];
                ['variance "size": \mu_2^{iso}/MD^2 = ' num2str(dps.vdison, 2)]
                ['root-mean-square "shape":  ',char(181),'FA = ' num2str(dps.ufa, 2)]
            };
        case 'Szczepankiewicz15'    
            title_str2 = {...
                ['mean "size": MD = ' num2str(dps.mdiso/1e-9, 2) '\cdot10^{-9} m^2/s'];
                ['mean-square "shape": \itV\rm_a/MD^2 = ' num2str(4/5*dps.msdanison, 2)];
                ['variance "size": \itV\rm_i/MD^2 = ' num2str(dps.vdison, 2)]
                ['root-mean-square "shape":  ',char(181),'FA = ' num2str(dps.ufa, 2)]
            };
        case 'Westin16 "QTI"'    
            title_str2 = {...
                ['mean "size": MD = ' num2str(dps.mdiso/1e-9, 2) '\cdot10^{-9} m^2/s'];
                ['mean-square "shape": C_\mu = ' num2str(dps.ufa^2, 2)];
                ['variance "size": C_{MD} = ' num2str(dps.vdison, 2)]
                ['root-mean-square "shape":  ',char(181),'FA = ' num2str(dps.ufa, 2)]
            };
        case 'Szczepankiewicz16 "DIVIDE"'    
            title_str2 = {...
                ['mean "size": MD = ' num2str(dps.mdiso/1e-9, 2) '\cdot10^{-9} m^2/s'];
                ['mean-square "shape": MK_a = ' num2str(dps.MKa, 2)];
                ['variance "size": MK_i = ' num2str(dps.MKi, 2)]
                ['root-mean-square "shape":  ',char(181),'FA = ' num2str(dps.ufa, 2)]
            };
        case 'Topgaard19'    
            title_str2 = {...
                ['mean "size": E[\itD\rm_{iso}] = ' num2str(dps.mdiso/1e-9, 2) '\cdot10^{-9} m^2/s'];
                ['mean-square "shape": E[(\itD\rm_A-\itD\rm_R)^2]/E[\itD\rm_{iso}]^2 = ' num2str(dps.msdanison, 2)];
                ['variance "size": Var[\itD\rm_{iso}]/E[\itD\rm_{iso}]^2 = ' num2str(dps.vdison, 2)]
            };
        otherwise
            title_str2 = {};
    end

    title(axh, title_str2,'FontSize',opt.mplot.fs,'FontWeight','normal')

end
