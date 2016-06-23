function mgui_analysis_plot(model_name, S, xps, h_top, h_bottom)


plot_fun_name = [model_name '_plot'];
plot_fn = [plot_fun_name '.m'];

if (~exist(plot_fn, 'file'))
    
    plot_err_msg([strrep(plot_fn, '_', '\_') ' not found']); 
    return;
    
end

% for convenient format in the models
S = mean(S, 1)';

try
    
    feval(plot_fun_name, S, xps, h_top, h_bottom);
    
catch me
    
    try
        
        feval(plot_fun_name, S, xps, h_top);
        
    catch me2
        
        plot_err_msg(me2.message);
        return;
        
    end
    
end


    function plot_err_msg(str)
        text(0,0,str, 'parent', h_top);
        ylim(h_top, [-1 1] * 10);
        axis(h_top, 'off');
    end


end