function [I_B,sc] = mio_rescale(I_B, I_A, M)
% function [I_B,sc] = mio_rescale(I_B, I_A, M)
%
% rescale B to A


sc = zeros(1,size(I_A,4));

for c = 1:numel(sc)
    
    tmp_a = I_A(:,:,:,c);
    tmp_b = I_B(:,:,:,c);
    
    if (all(tmp_a(:) == 0))
        error('undefined behaviour');
    end
    
    sc(c) = ...
        double(tmp_a(M(:) & ~isnan(tmp_b(:)))) \ ...
        double(tmp_b(M(:) & ~isnan(tmp_b(:))));
    
    I_B(:,:,:,c) = tmp_b / abs(sc(c));
end
