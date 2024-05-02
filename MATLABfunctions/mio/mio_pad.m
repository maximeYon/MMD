function I = mio_pad(I, pad_xyz)
% function I = mio_pad(I, pad_xyz)
%
% Adds zero around 'I' according to the three numbers in 'pad_xyz', 
% which may be negative (this leads to trimming instead)

if (all(pad_xyz == 0)), return; end

s = size(I);

for c = 1:3
    I = pad(I, pad_xyz(c), c);
end

    function I = pad(I, n, dim)
        
        if (n == 0), return; end
        
        if (n > 0)
            s = size(I);
            s(dim) = n;
            Z = zeros(s, class(I));
            I = cat(dim, Z, I, Z);
        else
            n = -n;
            
            ind = (1 + n):(size(I,dim)-n);
            
            switch (dim)
                case 1
                    I = I(ind,:,:,:);
                case 2
                    I = I(:,ind,:,:);
                case 3
                    I = I(:,:,ind,:);
            end
        end
    end
end