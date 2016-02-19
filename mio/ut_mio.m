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
        
end
