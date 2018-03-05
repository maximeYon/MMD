function [R_out, wgwf] = cfa_find_worst_gwf_rot(gwf, rf, dt, ips, mode, opt)

if nargin < 6
    opt = optimset('fminsearch');
    opt.MaxFunEvals = 5000;
    opt.MaxIter     = 5000;
    opt.TolX = 1e-3;
    opt.TolFun = 1e-3;
end

x0 = rand(1,3)*pi;

x = fminsearch(@(x0) wf_eval_func(x0), x0, opt);

R_out = cfa_euler_ang_to_rotmat(x(1), x(2), x(3));

wgwf  = gwf*R_out;

    function c = wf_eval_func(v)
        
        R = cfa_euler_ang_to_rotmat(v(1), v(2), v(3));
        
        [cc, cs, cp] = cfa_maxwell_bias(gwf*R, rf, dt, ips);
        
        switch mode
            case 1
                c = min(cc(:));
            case 2
                c = min(cs(:));
            case 3
                c = min(cp(:));
        end
        
    end

end




