function p = mio_volume_loop(fun, I, M, opt)
% function p = mio_volume_loop(fun, I, M, opt)
%
% Applies the function fun to all voxels where M is positive
%
% In order to execute the loop in parallel, initiate parpool before
% executing this function. If number of workers > 1, volume loop will be
% executed in parallel

if (nargin < 3), error('all inputs required'); end
if (nargin < 4), opt.present = 1; end

opt = mio_opt(opt);

% Init the output structure by fitting once to a voxel in the middle
%  (note that the output will be permuted before returning)
n_param = numel(fun(double(squeeze(...
    I(round(end/2), round(end/2), round(end/2), :)))));
p = zeros(n_param, size(I,1), size(I,2), size(I,3));

try % start parpool outside this function to have it run in parfor mode
    tmp = gcp('nocreate');
    n_workers = tmp.NumWorkers;
catch
    n_workers = 1;
end

% Loop over the data
fprintf('Starting the loop. Workers = %i\n', n_workers); 
for k = 1:size(I,3)
    fprintf('k=%3.0i', k);
    for j = 1:size(I,2)
        
        if (any(M(:,j,k)))            
            if (n_workers == 1) % single core
                for i = 1:size(I,1)
                    if (M(i,j,k) == 0), continue; end
                    p(:,i,j,k) = fun(double(squeeze(I(i,j,k,:))), 1);
                end
                
            else % do parallel
                parfor i = 1:size(I,1)
                    if (M(i,j,k) == 0), continue; end
                    p(:,i,j,k) = fun(double(squeeze(I(i,j,k,:))), 1);
                end
            end
        end

        % print a status report
        if (any(M(:,j,k))), marker = 'o'; else marker = '.'; end
        if (mod(j,4) == 0), fprintf(marker); end
    end
    disp(';');
end

end

p = permute(p, [2 3 4 1]);