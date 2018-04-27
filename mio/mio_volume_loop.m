function p = mio_volume_loop(fun, I, M, opt, S)
% function p = mio_volume_loop(fun, I, M, opt, S)
%
% Applies the function fun to all voxels where M is positive
%
% In order to execute the loop in parallel, initiate parpool before
% executing this function. If number of workers > 1, volume loop will be
% executed in parallel
%
% fun - accepts one or two input arguments (signal vector + optionally
%          supplementary data)
%
% I   - 4D matrix
%
% M   - 3D mask
%
% opt - option structure
%
% S   - 4D supplementary data
tic
if (nargin < 3), error('all inputs required'); end
if (nargin < 4), opt.present = 1; end
if (nargin < 5), S = []; end

opt = mio_opt(opt);

% Create functions w and w/o support information
if (~isempty(S))
    f_fun = @(x,y) fun(x,y);
    f_sup = @(K,i) K(:,i);
else
    f_fun = @(x,y) fun(x);
    f_sup = @(K,i) [];
end


% Expand 4D volume to 2D so that n x m is samples x voxels
siz = size(I);
I = reshape(I, prod(siz(1:3)), siz(4))';
M = reshape(M, prod(siz(1:3)),      1)';

% Retain only non-masked-out voxels.
si = find((M>0).*(~all(I==0,1))); clear M
I  = double(I(:,si));

if (~isempty(S))
    assert(all(siz==size(S)), 'S and I must be of equal size');
    S = reshape(S, prod(siz(1:3)), siz(4))';
    S = double(S(:,si));
end


% Init the output structure by fitting once to a voxel in the middle
%  (note that the output will be permuted before returning)
f = @(d) round(size(I,d)/2);
n_param = numel(f_fun( I(:, f(2)), f_sup(S, f(2)) ));


% Assume single thread, but check if parallel computing is wanted.
n_workers = 1;

if ~opt.mio.no_parfor
    try % start parpool outside this function to have it run in parfor mode
        tmp = gcp('nocreate');
        n_workers = tmp.NumWorkers;
    catch
        warning(['User requested parallel workers for mio_volume_loop but parpool is not active! '...
            'Continuing with job on a single local worker!'])
    end
end


if n_workers == 1
    
    out = zeros(n_param, size(I,2));
    
    for k = 1:size(I,2)
        out(:,k) = f_fun( I(:,k), f_sup(S,k) );
    end
    
else
    
    % Start parallel computing on all available workers using SPMD:
    % https://se.mathworks.com/help/distcomp/spmd.html
    % Also works for ONE worker, and ZERO workers execute the program
    % locally. SPMD was introduced in version Matlab V2008b. Wrapped
    % in "else" to support users without parallel computing toolbox.
    
    spmd (n_workers)
        I   = getLocalPart( codistributed(I) );
        S   = getLocalPart( codistributed(S) );
        out = zeros(n_param, size(I,2));  
        
        for k = 1:size(I,2)
            out(:,k) = f_fun( I(:,k), f_sup(S,k) );
        end
    end
    out = horzcat(out{:});
    
end
clear I

% Revert mask subsampling
p = zeros(n_param, prod(siz(1:3)));
p(:,si) = out;

% Revert 2D to 4D transform so that we get "x, y, z, parameter" dimensions
p = reshape(p', siz(1), siz(2), siz(3), n_param);

toc
