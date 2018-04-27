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

if (nargin < 3), error('all inputs required'); end
if (nargin < 4), opt.present = 1; end
if (nargin < 5), S = []; else assert(~isempty(S), 'S cannot be empty if provided'); end

opt = mio_opt(opt);

% Manage supplementary data
if (~isempty(S))
    f = @(d) size(I,d) == size(S, d);
    assert( f(1) & f(2) & f(3), 'S and I must be of equal size');
    
    % use this structure to manage the loop both with and without
    % supplementary information
    f_fun = @(x,y) fun(x,y);
    f_sup = @(i,j,k) squeeze(S(i,j,k,:));
else
    f_fun = @(x,y) fun(x);
    f_sup = @(i,j,k) [];
end

% Init the output structure by fitting once to a voxel in the middle
%  (note that the output will be permuted before returning)
f = @(d) round(size(I,d)/2);
n_param = numel(f_fun(double(squeeze(I(f(1), f(2), f(3), :))),...
    f_sup(f(1), f(2), f(3))));

% Expand 4D volume to 1+1D so that n x m is samples times voxels
siz = size(I);
I = reshape(I, prod(siz(1:3)), siz(4))';
M = reshape(M, prod(siz(1:3)),      1)';

% Retain only non-masked-out voxels.
si = find((M>0).*(~all(I==0,1))); clear M
I  = double(I(:,si));

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
        out(:,k) = f_fun( I(:,k) );
    end
    I = out;
    
else
    
    % Start parallel computing on all available workers
    % Here we use SPMD (https://se.mathworks.com/help/distcomp/spmd.html)
    % instead of parfor. Also works for a single worker. Zero workers execute
    % the program locally. SPMD was introduced in version 2008b. Wrapped in
    % "else" statemet to support users without parallel computing toolbox.
    
    spmd (n_workers)
        I = getLocalPart( codistributed(I) );
        out = zeros(n_param, size(I,2));
        
        for k = 1:size(I,2)
            out(:,k) = f_fun( I(:,k) );
        end
    end
    I = horzcat(out{:});
    
end
clear out

% Revert mask subsampling
p = zeros(n_param, prod(siz(1:3)));
p(:,si) = I;

% Revert 1+1D to 4D transform so that we get "x, y, z, parameter" dimensions
p = reshape(p', siz(1), siz(2), siz(3), n_param);
