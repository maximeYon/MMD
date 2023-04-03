clear all

stats_fn = '/Users/daniel/Dropbox/NMRdata/United/Yuan20200320_stats.xlsx';

stats_table = readtable(stats_fn);
sz = size(stats_table);

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

pred_array = zeros(sz(1),Nvars);
for nvar = 1:Nvars
    variable = variables{nvar};
    values = stats_table.(variable);
    
    pvalue = anova1(values,groups,'off');
%     display([variable ' ' num2str(pvalue)])
    
    anova_struct.variable{nvar} = variable;
    anova_struct.pvalue{nvar} = pvalue;
    
    pred_array(:,nvar) = values;
end

anova_table = struct2table(anova_struct);
anova_table = sortrows(anova_table,'pvalue','ascend');

resp = strcmp(groups,'HGG');
ind = find(isfinite(values));

resp = resp(ind);
pred_array = pred_array(ind,:);

%%

combinations = nchoosek(1:Nvars,2);
Ncombs = size(combinations,1);

emptycell = cell(Ncombs,1);
roc_struct.variable1 = emptycell;
roc_struct.variable2 = emptycell;
roc_struct.auclog = zeros(Ncombs,1);
roc_struct.aucsvm = zeros(Ncombs,1);
roc_struct.aucnb = zeros(Ncombs,1);

for ncomb = 1:Ncombs
% for ncomb = 1:100
    ind_var = combinations(ncomb,:);
    pred = pred_array(:,ind_var);

    mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
    score_log = mdl.Fitted.Probability; % Probability estimates


    % Compute the standard ROC curve using the probabilities for scores.

    [Xlog,Ylog,Tlog,AUClog] = perfcurve(resp,score_log,'true');

    % Train an SVM classifier on the same sample data. Standardize the data.

    mdlSVM = fitcsvm(pred,resp,'Standardize',true);

    % Compute the posterior probabilities (scores).

    mdlSVM = fitPosterior(mdlSVM);
    [~,score_svm] = resubPredict(mdlSVM);

    % The second column of score_svm contains the posterior probabilities of bad radar returns.


    % Compute the standard ROC curve using the scores from the SVM model.

    [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(resp,score_svm(:,mdlSVM.ClassNames),'true');

    % Fit a naive Bayes classifier on the same sample data.

    mdlNB = fitcnb(pred,resp);

    % Compute the posterior probabilities (scores).

    [~,score_nb] = resubPredict(mdlNB);

    % Compute the standard ROC curve using the scores from the naive Bayes classification.

    [Xnb,Ynb,Tnb,AUCnb] = perfcurve(resp,score_nb(:,mdlNB.ClassNames),'true');

    % Plot the ROC curves on the same graph.

%     figure
%     plot(Xlog,Ylog)
%     hold on
%     plot(Xsvm,Ysvm)
%     plot(Xnb,Ynb)
%     legend('Logistic Regression','Support Vector Machines','Naive Bayes','Location','Best')
%     xlabel('False positive rate'); ylabel('True positive rate');
%     title('ROC Curves for Logistic Regression, SVM, and Naive Bayes Classification')
%     hold off

    % Although SVM produces better ROC values for higher thresholds, logistic regression is usually better at distinguishing the bad radar returns from the good ones. The % ROC curve for naive Bayes is generally lower than the other two ROC curves, which indicates worse in-sample performance than the other two classifier methods.

    % Compare the area under the curve for all three classifiers.

    roc_struct.variable1{ncomb} = variables{ind_var(1)};
    roc_struct.variable2{ncomb} = variables{ind_var(2)};
    roc_struct.auclog(ncomb) = AUClog;
    roc_struct.aucsvm(ncomb) = AUCsvm;
    roc_struct.aucnb(ncomb) = AUCnb;
%     display([variables{ind_var(1)} ' ' variables{ind_var(2)} ' ' num2str(AUClog) ' ' num2str(AUCsvm) ' ' num2str(AUCnb)])

end
%%
roc_table = struct2table(roc_struct);
roc_table = sortrows(roc_table,'auclog','descend');
roc_table(1:100,:)
roc_table = sortrows(roc_table,'aucsvm','descend');
roc_table(1:100,:)
roc_table_2var = roc_table;
%%
% return
% for nvars = 1:Nvars
%     variables_temp = anova_table.variable(1:nvars);
%     ind_var = zeros(size(anova_table.variable));
%     for nvar = 1:nvars
%         ind_var = any([ind_var strcmp(variables',variables_temp{nvar})],2);
%     end
%     pred = pred_array(:,ind_var);
% 
%     mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
%     score_log = mdl.Fitted.Probability; % Probability estimates
% 
% 
%     % Compute the standard ROC curve using the probabilities for scores.
% 
%     [Xlog,Ylog,Tlog,AUClog] = perfcurve(resp,score_log,'true');
% 
%     % Train an SVM classifier on the same sample data. Standardize the data.
% 
%     mdlSVM = fitcsvm(pred,resp,'Standardize',true);
% 
%     % Compute the posterior probabilities (scores).
% 
%     mdlSVM = fitPosterior(mdlSVM);
%     [~,score_svm] = resubPredict(mdlSVM);
% 
%     % The second column of score_svm contains the posterior probabilities of bad radar returns.
% 
% 
%     % Compute the standard ROC curve using the scores from the SVM model.
% 
%     [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(resp,score_svm(:,mdlSVM.ClassNames),'true');
% 
%     % Fit a naive Bayes classifier on the same sample data.
% 
%     mdlNB = fitcnb(pred,resp);
% 
%     % Compute the posterior probabilities (scores).
% 
%     [~,score_nb] = resubPredict(mdlNB);
% 
%     % Compute the standard ROC curve using the scores from the naive Bayes classification.
% 
%     [Xnb,Ynb,Tnb,AUCnb] = perfcurve(resp,score_nb(:,mdlNB.ClassNames),'true');
% 
%     % Plot the ROC curves on the same graph.
% 
%     figure
%     plot(Xlog,Ylog)
%     hold on
%     plot(Xsvm,Ysvm)
%     plot(Xnb,Ynb)
%     legend('Logistic Regression','Support Vector Machines','Naive Bayes','Location','Best')
%     xlabel('False positive rate'); ylabel('True positive rate');
%     title('ROC Curves for Logistic Regression, SVM, and Naive Bayes Classification')
%     hold off
% 
%     % Although SVM produces better ROC values for higher thresholds, logistic regression is usually better at distinguishing the bad radar returns from the good ones. The % ROC curve for naive Bayes is generally lower than the other two ROC curves, which indicates worse in-sample performance than the other two classifier methods.
% 
%     % Compare the area under the curve for all three classifiers.
% 
%     AUClog
%     AUCsvm
%     AUCnb  
%     pause
% end
% 
% return
% %%
% ind_var = find(any([strcmp(variables,'dtd_gamma_MKt_90')
%    strcmp(variables,'dtd_gamma_MKi_75')
%    strcmp(variables,'dtd_fractions_2_10')
%    strcmp(variables,'dtd_msddelta_bin5_gray_75')
%    strcmp(variables,'dtd_vsddelta_90')
%    strcmp(variables,'dtd_gamma_MKa_75')
%    ],1));
% % ind_var = find(any([strcmp(variables,'dtd_gamma_MKt_50')],1));
% % ind_var = find(any([strcmp(variables,'dtd_gamma_MKi_50')
% %     strcmp(variables,'dtd_gamma_MKa_50')],1));
% % 
%     pred = pred_array(:,ind_var);
% 
%     mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
%     score_log = mdl.Fitted.Probability; % Probability estimates
% 
% 
%     % Compute the standard ROC curve using the probabilities for scores.
% 
%     [Xlog,Ylog,Tlog,AUClog] = perfcurve(resp,score_log,'true');
% 
%     % Train an SVM classifier on the same sample data. Standardize the data.
% 
%     mdlSVM = fitcsvm(pred,resp,'Standardize',true);
% 
%     % Compute the posterior probabilities (scores).
% 
%     mdlSVM = fitPosterior(mdlSVM);
%     [~,score_svm] = resubPredict(mdlSVM);
% 
%     % The second column of score_svm contains the posterior probabilities of bad radar returns.
% 
% 
%     % Compute the standard ROC curve using the scores from the SVM model.
% 
%     [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(resp,score_svm(:,mdlSVM.ClassNames),'true');
% 
%     % Fit a naive Bayes classifier on the same sample data.
% 
%     mdlNB = fitcnb(pred,resp);
% 
%     % Compute the posterior probabilities (scores).
% 
%     [~,score_nb] = resubPredict(mdlNB);
% 
%     % Compute the standard ROC curve using the scores from the naive Bayes classification.
% 
%     [Xnb,Ynb,Tnb,AUCnb] = perfcurve(resp,score_nb(:,mdlNB.ClassNames),'true');
% 
%     % Plot the ROC curves on the same graph.
% 
%     figure
%     plot(Xlog,Ylog)
%     hold on
%     plot(Xsvm,Ysvm)
%     plot(Xnb,Ynb)
%     legend('Logistic Regression','Support Vector Machines','Naive Bayes','Location','Best')
%     xlabel('False positive rate'); ylabel('True positive rate');
%     title('ROC Curves for Logistic Regression, SVM, and Naive Bayes Classification')
%     hold off
% 
%     % Although SVM produces better ROC values for higher thresholds, logistic regression is usually better at distinguishing the bad radar returns from the good ones. The % ROC curve for naive Bayes is generally lower than the other two ROC curves, which indicates worse in-sample performance than the other two classifier methods.
% 
%     % Compare the area under the curve for all three classifiers.
% 
%     AUClog
%     AUCsvm
%     AUCnb
% 
% return
emptycell = cell(Nvars,1);
clear roc_struct
roc_struct.variable = emptycell;
roc_struct.auclog = zeros(Nvars,1);
roc_struct.aucsvm = zeros(Nvars,1);
roc_struct.aucnb = zeros(Nvars,1);

for nvar = 1:Nvars
    variable = variables{nvar};
    values = stats_table.(variable);

%     pred = pred_array(:,nvar);
    pred = values(ind);

    mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
    score_log = mdl.Fitted.Probability; % Probability estimates


    % Compute the standard ROC curve using the probabilities for scores.

    [Xlog,Ylog,Tlog,AUClog] = perfcurve(resp,score_log,'true');

    % Train an SVM classifier on the same sample data. Standardize the data.

    mdlSVM = fitcsvm(pred,resp,'Standardize',true);

    % Compute the posterior probabilities (scores).

    mdlSVM = fitPosterior(mdlSVM);
    [~,score_svm] = resubPredict(mdlSVM);

    % The second column of score_svm contains the posterior probabilities of bad radar returns.


    % Compute the standard ROC curve using the scores from the SVM model.

    [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(resp,score_svm(:,mdlSVM.ClassNames),'true');

    % Fit a naive Bayes classifier on the same sample data.

    mdlNB = fitcnb(pred,resp);

    % Compute the posterior probabilities (scores).

    [~,score_nb] = resubPredict(mdlNB);

    % Compute the standard ROC curve using the scores from the naive Bayes classification.

    [Xnb,Ynb,Tnb,AUCnb] = perfcurve(resp,score_nb(:,mdlNB.ClassNames),'true');

    % Plot the ROC curves on the same graph.

%     figure
%     plot(Xlog,Ylog)
%     hold on
%     plot(Xsvm,Ysvm)
%     plot(Xnb,Ynb)
%     legend('Logistic Regression','Support Vector Machines','Naive Bayes','Location','Best')
%     xlabel('False positive rate'); ylabel('True positive rate');
%     title('ROC Curves for Logistic Regression, SVM, and Naive Bayes Classification')
%     hold off
% 
%     % Although SVM produces better ROC values for higher thresholds, logistic regression is usually better at distinguishing the bad radar returns from the good ones. The % ROC curve for naive Bayes is generally lower than the other two ROC curves, which indicates worse in-sample performance than the other two classifier methods.
% 
%     % Compare the area under the curve for all three classifiers.
% 
%     AUClog
%     AUCsvm
%     AUCnb
    roc_struct.variable{nvar} = variable;
    roc_struct.auclog(nvar) = AUClog;
    roc_struct.aucsvm(nvar) = AUCsvm;
    roc_struct.aucnb(nvar) = AUCnb;
%     display([variable ' ' num2str(AUClog) ' ' num2str(AUCsvm) ' ' num2str(AUCnb)])


end

roc_table = struct2table(roc_struct);
roc_table = sortrows(roc_table,'aucsvm','descend');

roc_table_1var = roc_table;
