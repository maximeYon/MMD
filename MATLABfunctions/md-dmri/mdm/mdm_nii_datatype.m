function datatype = mdm_nii_datatype(value, t)
% function datatype = mdm_nii_datatype(value, t)
%
% Returns the datatype of a nifti based on the header info in a 
% matlab-friendly format
%
% set t = 1 to change float64 --> double and float32 --> single

if (nargin < 2), t = 0; end 

switch (value)
    case 2048
        datatype = 'complex256';
        error('datatype not supported');
    case 1792
        datatype = 'complex128';
    case 1536
        datatype = 'float128';
        error('datatype not supported');
    case 1280
        datatype = 'uint64';
    case 1024
        datatype = 'int64';
    case 512
        datatype = 'uint16';
    case 256
        datatype = 'int8';
    case 128
        error('datatype rgb, should be caught by 24 bits per voxel')
    case 64
        if (t)
            datatype = 'double';
        else
            datatype = 'float64';
        end
    case 32
        datatype = 'complex64';
    case 16
        if (t)
            datatype = 'single';
        else
            datatype = 'float32';
        end
    case 8
        datatype = 'int32';
    case 4
        datatype = 'int16';
    case 2
        datatype = 'uint8';
    case 1
        if (t)
            error('bit1 not defined in matlab');
        else
            datatype = 'bit1'; % one bit
        end
    otherwise
        error('could not determine datatype correctly');
end

