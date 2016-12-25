function s = mdm_powder_average(varargin)
% function s = mdm_powder_average(s, o_path, opt)
%
% Average over rotations. Image volumes with identical rotations is defined
% from s.xps.a_ind
%
% To do: find a way of keeping track of number of averages per step

warning('legacy, use mdm_s_powder_average instead');

s = mdm_powder_average(varargin{:});