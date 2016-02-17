function mfs = mio_volume_loop(fun, I, M, mfs)
% function mfs = mio_fit_model(s, fun, opt)

if (nargin < 3), error('all inputs required'); end
if (nargin < 4), mfs = []; end % fill in an existing mfs on own risk

% Create the model fit structure
m = fun([], 0);

for f = fieldnames(m)'
    mfs.(f{1}) = zeros([size(I,1) size(I,2) size(I,3) numel(m.(f{1}))]);
end

% Loop over the data
disp('Starting the loop');
for k = 1:size(I,3)
    fprintf('k=%i', k);
    for j = 1:size(I,2)
        
        marker = '.';
        for i = 1:size(I,1)
            
            if (M(i,j,k) == 0), continue; end
            marker = 'o';
            
            m = fun(squeeze(I(i,j,k,:)), 1);
            
            for f = fieldnames(m)'
                mfs.(f{1})(i,j,k,:) = m.(f{1});
            end
        end        
        if (mod(j,4) == 0), fprintf(marker); end
    end
    disp(';');
end

mfs.M = M;