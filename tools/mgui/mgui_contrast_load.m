function [I, header, filename, ext, xps] = mgui_contrast_load(EG)
% function [I, header, filename, ext, xps] = mgui_contrast_load(EG)
%
% Return the image volume and filename based on the select struct

do_flip = 1;

filename = EG.browse.filename;

[~,~,ext] = msf_fileparts(filename);

% Search for .gz files if the orignal file is not present
if (~exist(filename, 'file') && exist([filename '.gz'], 'file'))
    filename = [filename '.gz'];
end

% Nifti files are in LPS (some, at least, check this!)
% Matlab images are displayed as YX
[I, header] = mdm_nii_read(filename);
I = double(I);

header.my_hdr.ori = mdm_nii_oricode(header);

[~,name] = fileparts(filename);
header.my_hdr.protocol_name = name(max(find(name == '_', 2))+1:end);

if (header.bitpix == 24)
    I = permute(I, [2 3 4 1]);
end

if (header.bitpix == 8) && (header.dim(2) == 3)
    I = permute(I, [1 2 3 4]);
    I = reshape(I(:), 128, 128, 60, 3);
    header.bitpix = 24;
end

if (header.scl_slope ~= 0)
    I = I * double(header.scl_slope);
end

if (header.scl_inter ~= 0)
    I = I + double(header.scl_inter);
end

% Possibly scale according to my_header.ref_scale
if (isfield(header.my_hdr, 'ref_scale') && ...
        header.my_hdr.ref_scale ~= 0)
    I = I / header.my_hdr.ref_scale;
end

% Check if there is an xps-file available
xps_filename = mdm_xps_fn_from_nii_fn(filename);
if (exist(xps_filename, 'file'))
    xps = mdm_xps_load(xps_filename);
else
    xps = [];
end

% If there's no xps, check for gdir
if (isempty(xps))
    gdir_filename = mdm_fn_nii2gdir(filename);
    
    if (exist(gdir_filename, 'file'))
        xps = mdm_xps_from_gdir(gdir_filename);
    end
end

    
% Rotate the volume
if (do_flip) && (isfield(header, 'my_hdr') && (isfield(header.my_hdr,'ori')))
    [I,p] = mgui_misc_flip_volume(I, header.my_hdr.ori, EG.conf.ori);
    
    if (isfield(header, 'pixdim'))
        header.pixdim(2:4) = header.pixdim(p + 1);
    end
end


end
