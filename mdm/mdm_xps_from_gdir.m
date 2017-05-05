function xps = mdm_xps_from_gdir(gdir_fn, delimeter, b_delta)
% function xps = mdm_xps_from_gdir(gdir_fn, delimeter, b_delta)
%
% read a gradient textfile in the Lund format, which is defined as follows
%   n_x, n_y, n_z, b
% 
% units of b are in s/mm2, i.e., b = 1000 s/mm2 for a standard DTI
%
% the xps is always in SI units
%
% optional arguments
% delimeter, asumed to be ',' unless specified
% b_delta, assumed to be 1 unless specified (e.g. PGSE / LTE / SDE)

if (nargin < 2), delimeter = []; end
if (nargin < 3), b_delta = 1; end

if (isempty(delimeter)), delimeter = ','; end
if (b_delta < -0.5) || (b_delta > 1), error('b_delta out of range'); end

if (~exist(gdir_fn)), error('Did not find file %s', gdir_fn); end

% read text and split by delimeter, keeping only non-empty entries to solve
% problems of multiple consequtive delimeters such as spaces
txt = mdm_txt_read(gdir_fn);
f = @(x) x(~cellfun(@isempty, x));

tmp = cellfun(@(x) strsplit(x, delimeter), txt, 'uniformoutput', 0);
tmp = cellfun(f, tmp, 'uniformoutput', 0);
tmp = cellfun(@(x) str2double(x)', tmp, 'uniformoutput', 0);
tmp = cell2mat(tmp);

% compile the output
xps.b   = tmp(4,:)' * 1e6; 
xps.u   = tmp(1:3,:)';

% demand normalization of the vectors
xps.u   = xps.u ./ repmat(eps + sqrt(sum(xps.u.^2,2)), 1, 3);

xps.n   = numel(xps.b);
% xps.bt  = tm_1x3_to_1x6(xps.b, zeros(size(xps.b)), xps.u);
xps.bt  = tm_tpars_to_1x6(xps.b, b_delta, xps.u);
xps.bt2 = tm_1x6_to_1x21(xps.bt);
xps.b_delta = zeros(size(xps.b)) + b_delta;
xps.b_eta   = zeros(size(xps.b));

% 
% if (b_delta == 0)
%     xps.bt = repmat([1/3 1/3 1/3 0 0 0], size(xps.bt,1), 1) .* repmat(xps.b, 1, 6);
% elseif (b_delta ~= 1)
%     
%     for c = 1:size(xps.bt, 1)
%         
%         p = tm_3x3_to_tpars( tm_1x6_to_3x3( xps.bt(c, :) ) );
%         b = p.trace;
%     
%         b_stick = xps.bt(c,:) / b;
%         b_sphere = tm_3x3_to_1x6( eye(3,3) / 3 );
%     
%         % b = c1 + c2
%         c2 = b * b_delta;
%         c1 = b - c2;
%     
%         xps.bt(c,:) = c1 * b_sphere + c2 * b_stick;
%     
%     end
%     
% %     error('not yet implemented');
% end


