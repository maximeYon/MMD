function [n,par,perp,theta,phi,d0,rpar,rperp,r1,r2] = dtor1r2d_nodes2par(dtor1r2d_nodes)

if ~isempty(dtor1r2d_nodes)
    n = dtor1r2d_nodes(1);
    if n > 0
        m = numel(dtor1r2d_nodes(2:end))/n;
        dtor1r2d_nodes_array = reshape(dtor1r2d_nodes(2:end),[m n]);
        par = dtor1r2d_nodes_array(1,:)';
        perp = dtor1r2d_nodes_array(2,:)';
        theta = dtor1r2d_nodes_array(3,:)';
        phi = dtor1r2d_nodes_array(4,:)';
        d0 = dtor1r2d_nodes_array(5,:)';
        rpar = dtor1r2d_nodes_array(6,:)';
        rperp = dtor1r2d_nodes_array(7,:)';
        r1 = dtor1r2d_nodes_array(8,:)';
        r2 = dtor1r2d_nodes_array(9,:)';
    else
        n = 0;
        par = [];
        perp = [];
        theta = [];
        phi = [];
        d0 = [];        
        rpar = [];        
        rperp = []; 
        r1 = [];
        r2 = [];
    end
else
    n = 0;
    par = [];
    perp = [];
    theta = [];
    phi = [];
    d0 = [];
    rpar = [];
    rperp = [];
    r1 = [];
    r2 = [];
end
