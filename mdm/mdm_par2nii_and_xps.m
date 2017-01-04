function [nii_fn, xps] = mdm_par2nii_and_xps(par_fn, nii_fn)
% function [nii_fn, xps] = mdm_par2nii_and_xps(par_fn, nii_fn)


txt = mdm_txt_read(par_fn);
txt = txt(cellfun(@(x) ~isempty(strtrim(x)), txt)); % strip empty lines

% Read the header of the .par file
for c = find(cellfun(@(x) x(1) == '.', txt))
    
    % Find the header information in the par-file
    tline = txt{c};
    
    t = strsplit(tline(2:end),':');
    
    key = char(strtrim(t{1}));
    val = char(strtrim(t{2}));
    
    switch (key)
        case 'Patient name'
            s.subject_name = val;
        case 'Protocol name'
            s.protocol_name = val;
        case 'Repetition time [ms]'
            s.TR = str2num(val);
        case 'Scan resolution  (x, y)'
            s.image_size = str2num(val);
        case 'Max. number of diffusion values'
            s.n_samples = str2num(val); % warning: this assumption may not hold
        case 'sMax. number of gradient orients'
            s.n_rows = str2num(val);
        case 'Angulation midslice(ap,fh,rl)[degr]'
            s.angulation_midslice = str2num(val) * pi/180; % radians
        case 'Off Centre midslice(ap,fh,rl) [mm]'
            s.off_centre_midslice = str2num(val);
        case 'Preparation direction'
            s.preparation_direction = val;
        case 'FOV (ap,fh,rl) [mm]'
            s.fov = str2num(val);
    end
end

% Read the information specific for each image
txt2 = txt(cellfun(@(x) (x(1) ~= '.') && (x(1) ~= '.') && (x(1) ~= '#'), txt));
if (numel(txt2) == 0), error('strange PAR file'); end

for c = 1:numel(txt2)
    
    try
        info = eval(['[' txt2{c} ']']);
    catch me
        disp(me.message);
        error('cannot read line in PAR file');
    end
    
    s.info(c).c_slice           = info(1);
    s.info(c).c_echo            = info(2);
    s.info(c).c_dynamic         = info(3);
    s.info(c).c_cardiac         = info(4);
    s.info(c).image_type        = info(5);
    s.info(c).scanning_sequence = info(6);
    s.info(c).index_in_rec      = info(7);
    s.info(c).n_pixel_bits      = info(8);
    s.info(c).scan_percentage   = info(9);
    s.info(c).image_size        = info([10 11]);
    s.info(c).rescale_intercept = info(12);
    s.info(c).rescale_slope     = info(13);
    s.info(c).scale_slope       = info(14);
    s.info(c).window_center     = info(15);
    s.info(c).window_width      = info(16);
    s.info(c).angulation        = info([17 18 19]);
    s.info(c).offcentre         = info([20 21 22]);
    s.info(c).slice_thickness   = info(23);
    s.info(c).slice_gap         = info(24);
    s.info(c).display_orientation  = info(25);
    s.info(c).slice_orientation    = info(26);
    s.info(c).fmri_status          = info(27);
    s.info(c).image_type           = info(28);
    s.info(c).pixel_spacing        = info([29 30]);
    s.info(c).echo_time            = info(31);
    s.info(c).dyn_scan_begin       = info(32);
    s.info(c).trigger_time         = info(33);
    s.info(c).diffusion_b_factor   = info(34);
    s.info(c).n_averages           = info(35);
    s.info(c).image_flip_angle     = info(36);
    s.info(c).cardiac_bpm          = info(37);
    s.info(c).min_rr_interval      = info(38);
    s.info(c).max_rr_interval      = info(39);
    s.info(c).turbo_factor         = info(40);
    s.info(c).inversion_delay      = info(41);
    s.info(c).b_value_number       = info(42);
    s.info(c).gradient_orientation = info(43);
    s.info(c).contrast_type        = info(44);
    s.info(c).diffusion_anisotropy = info(45);
    s.info(c).diffusion_dir        = info([46 47 48]);
end

if (~(sum(std([s.info(:).image_size],1,2)) == 0))
    error('This routine is not made for images of different sizes in par/rec');
end

% Calculate some information
s.n_images = c;
s.n_slice = max([s.info(:).c_slice]);
s.n_dynamic = max([s.info(:).c_dynamic]);

% Read REC data
rec_fn = [par_fn(1:find(par_fn == '.', 1, 'first')) 'REC'];
if (~exist(rec_fn, 'file'))
    error(['Data (.rec) file does not exist (' rec_fn ')']);
end

fid = fopen(rec_fn,'r');
I = fread(fid, inf, 'ushort');
fclose(fid);

n_x = s.info(end).image_size(1);
n_y = s.info(end).image_size(2);
n_z = s.n_slice;

I = reshape(I, n_x, n_y, numel(I) / (n_x * n_y));

% Weed out geomean
good_ind = ~ ( ([s.info.diffusion_b_factor] > 0) & ...
    (sqrt(sum(reshape([s.info.diffusion_dir],3,s.n_images).^2,1)) == 0) );

if (sum(good_ind) == 0), good_ind = 1:size(s.I,3); end % never remove all

I = I(:,:,good_ind);
s.info = s.info(good_ind); 

% Split data into separate images
n_vol = numel(I) / (n_x * n_y * n_z);
x = repmat(1:n_z, n_vol, 1);
if (all(x(:) == [s.info(:).c_slice]'))
    I = reshape(I, n_x, n_y, n_vol, n_z);
    I = permute(I, [1 2 4 3]);
else
    I = reshape(I, n_x, n_y, n_z, n_vol);
end

s.info = s.info([s.info.c_slice] == 1); % assume first slice is sufficient

% Rescale images
%
% # === PIXEL VALUES =============================================================
% #  PV = pixel value in REC file, FP = floating point value, DV = displayed value on console
% #  RS = rescale slope,           RI = rescale intercept,    SS = scale slope
% #  DV = PV * RS + RI             FP = DV / (RS * SS)
for c = 1:size(I,4)
    PV = double(I(:,:,:,c));
    DV = PV * s.info(c).rescale_slope + s.info(c).rescale_intercept;
    FP = DV / (s.info(c).rescale_slope * s.info(c).scale_slope);
    I(:,:,:,c) = FP;
end





% Construct nifti header
h = mdm_nii_h_empty;

h.pixdim(2:4) = [s.info(1).pixel_spacing s.info(1).slice_thickness];
h.dim(2:4)    = [n_x n_y n_z];

warning('assuming LAS orientation for nifti, which may be wrong');
I = I(:,end:-1:1,:,:); % flip in to LAS
h.sform_code    = 1;
h.srow_x        = [-h.pixdim(2) 0 0 +single(h.dim(2)) * h.pixdim(2) / 2];
h.srow_y        = [0 +h.pixdim(3) 0 -single(h.dim(3)) * h.pixdim(3) / 2];
h.srow_z        = [0 0 +h.pixdim(4) -single(h.dim(4)) * h.pixdim(4) / 2];

mdm_nii_write(I, nii_fn, h);




% Construct xps
xps.n = size(I,4);

xps.b       = [s.info.diffusion_b_factor]' * 1e6;
xps.u       = reshape([s.info.diffusion_dir],3,xps.n)';
xps.s_ind   = [s.info.c_dynamic]';
xps.b_ind   = [s.info.b_value_number]';
xps.br_ind  = [s.info.gradient_orientation]';
xps.te      = [s.info.echo_time]' * 1e-3;
xps.bt      = tm_1x3_to_1x6(xps.b, zeros(size(xps.b)), xps.u);



