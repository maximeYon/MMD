function xps = mdm_philips_fexi_prot2xps(prot_fn, xps)
% function xps = mdm_philips_fexi_prot2xps(prot_fn, xps)
%
% Parse a Philips protocol text file in prot_fn, and store it in an xps
%
% The protocol file must have been generated with the FEXI patch
%
% An xps must be provided and some fields must have been filled in 
% (see this file for details)


% associate tm so different s_ind's
if (~isfield(xps, 's_ind')), error('xps.s_ind required'); end 
if (~isfield(xps, 'mde_b2_ind')), error('xps.mde_b2_ind required'); end 
    

% Read protocol and check that it is a FEXI protocol
p = mdm_philips_parse_txt_protocol(prot_fn);
f = fieldnames(p);

if (~any(cellfun(@(x) ~isempty(strfind('fT2prep__DWFILTER', x)), f)))
    error('this does not seem to be a fexi protocol file (%s)', prot_fn);
end

if (~strcmpi(p.fT2prep__DWFILTER, 'yes'))
    error('The DW filter did was not active, not a FEXI acquisition');
end

% Extract useful information
n_tm = p.f____dyn_scans;        

if (n_tm == 1), error('only one mixing time, not a proper FEXI'); end
if (max(xps.s_ind) ~= n_tm), error('xps.s_ind should have a range from 1 to %i', n_tm); end


% filter active for s_ind > 1, i.e., all but the first dynamic
xps.mde_b1 = p.f____filter_b_value_smm2 * 1e6 * (xps.s_ind > 1);

% try to use existing b-values for the mde_b2
xps.mde_b2 = p.f____b_factors(xps.mde_b2_ind)' * 1e6;

% get mixing times
xps.mde_tm12 = p.f____filter_tm_ms(xps.s_ind)' * 1e-3;

% fill in for averaging
[~,~,xps.a_ind] = unique([xps.s_ind xps.mde_b2_ind], 'rows');


% Relevant fields        
% f____dyn_scans;
% f____filter_tm_ms
% f____dyn_scans: 4
% f____echo_time_ms: 39
% f__enable_filter: 'yes'
% f____slice_selective: 'yes'
% f____filter_tm_ms: [1x32 double]
% f____filter_delta_ms: 11
% f____filter_b_value_smm2: 840
% f____custom_direction: 'no'        
        
% DTI:
% fDiffusion_mode: 'DTI'
% f____sequence: 'SE'
% f____gradient_duration: 'maximum'
% f____gradient_overplus: 'yes'
% f____directional_resolution: 'low'
% f____nr_of_b_factors: 9
% f____b_factor_order: 'user defined'
% f____b_factors: [1x16 double]
% f____average_high_b: 'no'

