clear all

stats_fn = '/Users/daniel/Dropbox/NMRdata/United/Yuan20200320_stats.xlsx';

stats_table = readtable(stats_fn);
columns = stats_table.Properties.VariableNames;
grade_234 = stats_table.PathologyLabel;
grade_lh = cell(size(grade_234));
ind_l = strcmp(grade_234,'2');
ind_h = any([strcmp(grade_234,'3') strcmp(grade_234,'4')],2);
[grade_lh{ind_l}] = deal('LGG');
[grade_lh{ind_h}] = deal('HGG');

variables = columns([2 4:end]);
Nvars = numel(variables);
groups = grade_lh;

emptycell = cell(Nvars,1);
anova_struct.variable = emptycell;
anova_struct.pvalue = emptycell;

for nvar = 1:Nvars
    variable = variables{nvar};
    values = stats_table.(variable);
    
    pvalue = anova1(values,groups,'off');
%     display([variable ' ' num2str(pvalue)])
    
    anova_struct.variable{nvar} = variable;
    anova_struct.pvalue{nvar} = pvalue;
end

anova_table = struct2table(anova_struct);

anova_table = sortrows(anova_table,'pvalue','ascend');

% for nvar = 1:Nvars
%     variable = anova_table.variable{nvar};
%     values = stats_table.(variable);
%     
%     pvalue = anova1(values,groups,'on');
%     display([variable ' ' num2str(pvalue)])
%     pause(1)
%     delete(gcf), delete(gcf)
%     
% end
