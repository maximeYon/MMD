function str_out = mdm_bruker_str_pv2tex(str_in)

str_out = str_in;

ind_u = find(str_out=='_');
for n_u = 1:numel(ind_u)
    ind_end = numel(str_out);
    ind_pre = (1:(ind_u(n_u)-1)+n_u-1);
    ind_post = (ind_u(n_u)+n_u):ind_end;
    str_out = strcat(str_out(ind_pre),'\_',str_out(ind_post));
end
str_out(str_out=='<') = '';
str_out(str_out=='>') = '';


