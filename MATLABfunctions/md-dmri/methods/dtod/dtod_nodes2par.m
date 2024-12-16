function [n,par,perp,theta,phi,d0,rpar,rperp] = dtod_nodes2par(dtod_nodes)

if ~isempty(dtod_nodes)
    n = dtod_nodes(1);
    if n > 0
        m = numel(dtod_nodes(2:end))/n;
        dtod_nodes_array = reshape(dtod_nodes(2:end),[m n]);
        par = dtod_nodes_array(1,:)';
        perp = dtod_nodes_array(2,:)';
        theta = dtod_nodes_array(3,:)';
        phi = dtod_nodes_array(4,:)';
        d0 = dtod_nodes_array(5,:)';
        rpar = dtod_nodes_array(6,:)';
        rperp = dtod_nodes_array(7,:)';
    else
        n = 0;
        par = [];
        perp = [];
        theta = [];
        phi = [];
        d0 = [];        
        rpar = [];        
        rperp = [];        
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
end
