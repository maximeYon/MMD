function mgui_analysis_plot(model_name, S, xps, h_top, h_bottom)
% function mgui_analysis_plot(model_name, S, xps, h_top, h_bottom)


if (numel(S) == 0)
    mgui_analysis_plot_overview(S, xps, h_top, h_bottom);
    return;
end
    

plot_fun_name   = [model_name '_plot'];
fun_1d_data2fit = [model_name '_1d_data2fit'];
fun_1d_fit2data = [model_name '_1d_fit2data'];
fun_opt         = [model_name '_opt'];

% for convenient format in the models
MS = mean(S, 1)';

if (exist([plot_fun_name '.m'], 'file'))
    try % standard plot function

        feval(plot_fun_name, MS, xps, h_top, h_bottom);
        return;
        
    catch me
        fprintf('%s err: %s\n', plot_fun_name, me.message);
    end
    
    % if complex, try with magnitude data
    if (any(~isreal(MS(:))))
        try % standard plot function
            
            feval(plot_fun_name, abs(MS), xps, h_top, h_bottom);
            warning('Performed analysis on magnitude part of complex data ? dataset will perhaps not work with fit methods');
            return;
            
        catch me
            fprintf('%s err: %s\n', plot_fun_name, me.message);
        end
    end
    
    try % with one argument less
        
        feval(plot_fun_name, MS, xps, h_top);
        return;
        
    catch me
        fprintf('%s err: %s\n', plot_fun_name, me.message);
    end
end

try % standard 1d fit

    opt = feval(fun_opt);
    m = feval(fun_1d_data2fit, MS, xps, opt);
    S_fit = feval(fun_1d_fit2data, m, xps)';
    
    mgui_analysis_plot_overview(S, xps, h_top, h_bottom, S_fit);
    return;
catch me
    fprintf('%s err: %s\n', plot_fun_name, me.message);
end


str = 'Could not fit/plot selected model';

text(0,0,str, 'parent', h_top);
ylim(h_top, [-1 1] * 10);
axis(h_top, 'off');



