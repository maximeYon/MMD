function p = mio_volume_loop(fun, I, M, opt)
% function p = mio_volume_loop(fun, I, M, opt)

if (nargin < 3), error('all inputs required'); end
if (nargin < 4), opt.present = 1; end

opt = mio_opt(opt);

% Init the output structure (but permute before outputting)
p = zeros(numel(fun([], 0)), size(I,1), size(I,2), size(I,3));


% Run parallel?
do_delete_parpool = 0;
if (opt.mio.n_core > 1)
    if (isempty(gcp('nocreate')))
        parpool(opt.mio.n_core);
        do_delete_parpool = 1;
    end    
end

try
    tmp = gcp('nocreate');
catch 
    tmp = [];
end

if (isempty(tmp))
    n_workers = 1;
else
    n_workers = tmp.NumWorkers;
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
                    p(:,i,j,k) = fun(squeeze(I(i,j,k,:)), 1);
                end
                
            else % do parallel
                parfor i = 1:size(I,1)
                    if (M(i,j,k) == 0), continue; end
                    p(:,i,j,k) = fun(squeeze(I(i,j,k,:)), 1);
                end
            end
        end

        % print a status report
        if (any(M(:,j,k))), marker = 'o'; else marker = '.'; end
        if (mod(j,4) == 0), fprintf(marker); end
    end
    disp(';');
end

% close parpool if it was created
if (do_delete_parpool)
    delete(gcp('nocreate'));
end

p = permute(p, [2 3 4 1]);