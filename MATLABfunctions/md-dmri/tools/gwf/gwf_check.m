function gwf_check(gwf, rf, dt)
% function gwf_check(gwf, rf, dt)
%
% Check that a gwf is properly formatted

if (~ismatrix(gwf))
    error('ndims(gwf) > 2'); 
end

if (size(gwf,2) > 3), error('gwf should be n x (1, 2, or 3)'); end
if (size(rf, 2) ~= 1), error('rf should be n x 1'); end

if (size(gwf,1) ~= size(rf, 1))
    error('Expected size(gwf,1) == size(rf,1), got %i ~= %i', ...
        size(gwf, 1), size(rf, 1)); 
end


if (size(gwf,2) > 3)
    error('size(gwf) = %i x %i, should be n x (1,2,3)', ...
        size(gwf,1), size(gwf,2)); 
end

if (numel(rf) ~= size(gwf, 1)) || (size(rf,2) ~= 1)
    error('size(rf) = %i x %i, should be n x 1', ...
        size(rf,1), size(rf,2)); 
end    




