function [n,par,perp] = dtd_pa_nodes2par(dtd_nodes)

if ~isempty(dtd_nodes)
    n = dtd_nodes(1);
    m = numel(dtd_nodes(2:end))/n;
    dtd_nodes_array = reshape(dtd_nodes(2:end),[m n]);
    par = dtd_nodes_array(1,:)';
    perp = dtd_nodes_array(2,:)';
else
    n = 0;
    par = [];
    perp = [];
end
