function fn = ut_mio(c_ut)
% function fn = ut_mio(c_ut)
%
% Run unit tests on the files in this package

if (nargin == 0), fn = 1; return; end

switch (c_ut)
    
    case 1
        fn = 'mio_pad.m';
        
        I = zeros(10,10,10);
        
        I2 = mio_pad(I, [1 2 3]);
        
        if (~all(size(I2) == [12 14 16])), error('%s, ut_mio test %i', fn, c_ut); end
        
        I3 = mio_pad(I2, -[1 2 3]);
        
        if (~all(size(I3) == size(I))), error('%s, ut_mio test %i, step2', fn, c_ut); end
        
    case {2,3,4,5}
        fn = 'mio_volume_loop.m';
        
        opt.verbose = 0;
        
        I = ones(11,12,13,14);
        M = I(:,:,:,1) > 0;
        
        % allow test with and without parpool
        try % start parpool outside this function to have it run in parfor mode
            tmp = gcp('nocreate');
            n_workers = tmp.NumWorkers;
        catch
            n_workers = 1;
        end
        
        switch (c_ut)
            case 2
                opt.mio.no_parfor = 1;
                P = mio_volume_loop(@sum, I, M, opt);
                if (~all(P(:) == 14)), error('%s, ut_mio test %i', fn, c_ut); end
            case 3
                opt.mio.no_parfor = 1;
                P = mio_volume_loop(@(x,y) sum(x) + sum(y), I, M, opt, I);
                if (~all(P(:) == 14*2)), error('%s, ut_mio test %i', fn, c_ut); end
            case 4
                if (n_workers == 1), error('%s, ut_mio test %i; run parpool(x) first to test such functions', fn, c_ut); end
                P = mio_volume_loop(@sum, I, M, opt);
                if (~all(P(:) == 14)), error('%s, ut_mio test %i', fn, c_ut); end
                
            case 5
                if (n_workers == 1), error('%s, ut_mio test %i; run parpool(x) first to test such functions', fn, c_ut); end
                P = mio_volume_loop(@(x,y) sum(x) + sum(y), I, M, opt, I);
                if (~all(P(:) == 14*2)), error('%s, ut_mio test %i', fn, c_ut); end
                
        end
        
end
