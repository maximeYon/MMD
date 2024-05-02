function [Iout,ar,c_slice_eff] = mgui_misc_format_image(EG, I, c_slice, c_dim, hdr, c_volume, do_sum, do_average)
% function Iout = mgui_misc_format_image(EG, I, c_slice, c_dim, hdr)

% Init
if (nargin < 5), hdr = []; end
if (nargin < 6), c_volume = []; end
if (nargin < 7), do_sum = []; end
if (nargin < 8), do_average = []; end

if (isempty(hdr)), hdr.present = 1; end
if (isempty(c_volume)), c_volume = 1; end
if (isempty(do_sum)), do_sum = 0; end
if (isempty(do_average)), do_average = 0; end


Iout = 0;
ar = [1 1];

if (numel(I) == 0), return; end

if (ndims(I) == 4)
    
    % Is RBG uint8?
    if ((isfield(hdr,'bitpix') && (hdr.bitpix == 24) && (size(I,4) == 3)))
        
        c_slice_eff = max(1, min(max(c_slice), size(I, c_dim)));
        switch (c_dim)
            case 1
                I = (I(c_slice_eff,:,:,:));
            case 2
                I = (I(:,c_slice_eff,:,:));
            case 3
                I = (I(:,:,c_slice_eff,:));
        end

        [I2, ar] = mgui_misc_format_image(EG, I(:,:,:,1), 1, c_dim, hdr);
        I2 = repmat(I2, 1, 1, 3);
        I2(:,:,2) = mgui_misc_format_image(EG, I(:,:,:,2), 1, c_dim, hdr);
        I2(:,:,3) = mgui_misc_format_image(EG, I(:,:,:,3), 1, c_dim, hdr);
        
        
        if (EG.roi.I_max ~= 0)
            Iout = I2 / EG.roi.I_max * 1.2;
        end
        
        return;
    end
end


c_slice_eff = max(1, min(max(c_slice), size(I, c_dim)));

if (do_sum)
    
    % Possibly sum only over a slab
    if (numel(c_slice) > 1)
        switch (c_dim)
            case 1
                I = I(c_slice(1):c_slice(2),:,:,c_volume);
            case 2
                I = I(:,c_slice(1):c_slice(2),:,c_volume);
            case 3
                I = I(:,:,c_slice(1):c_slice(2),c_volume);
        end
    end
    
    Iout = squeeze(sum(I, c_dim));
    if (do_average && size(I, c_dim) > 0)
        Iout = Iout / size(I, c_dim);
    end
    
else
    
    switch (c_dim)
        case 1
            Iout = squeeze(I(c_slice_eff,:,:,c_volume));
        case 2
            Iout = squeeze(I(:,c_slice_eff,:,c_volume));
        case 3
            Iout = squeeze(I(:,:,c_slice_eff,c_volume));
    end
end

% Matlab image coordinate system
Iout = rot90(Iout);

% Calculate the aspect ratio
if (isstruct(hdr) && isfield(hdr, 'pixdim'))
    pixdim = hdr.pixdim(2:4);
    pixdim = pixdim(:)';
else
    pixdim = [ 1 1 1 ];
end
pixdim = pixdim(1:3 ~= c_dim);

ar = pixdim / min(pixdim);

if (any(isnan(ar))), ar = [1 1]; end
