function xps = jm_xps_from_log( diff_array )

if( size(diff_array,2) < 5 )
    error('ERROR: diff_array does not contain sufficient information.');
end

% b-Value for each measurement in SI
xps.b = diff_array(:, 1) * 1e+06;
% Total number of measurements
xps.n = size(diff_array, 1);
% Gradient directions and b-tensor
gdir = diff_array(:, 2:4);
lambda_axial = diff_array(:,5);
b_delta = lambda_axial - (ones(size(lambda_axial)) - lambda_axial)/2;

bt = tm_tpars_to_1x6( xps.b, b_delta, gdir );
xps = mdm_xps_from_bt(bt);

xps.u = gdir;
