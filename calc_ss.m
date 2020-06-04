function [acc,acc_m,acc_f,acc_hcp,acc_m_hcp,acc_f_hcp] = calc_ss(samplesize_th,FC_UKB,cov,FC1_regressed,FC2_regressed,FC3_regressed,FC4_regressed,Gender_hcp)
for rpt = 1:100
    clear FC_regressed
    
    % random selecting subjects in order to keep the balance between female and male subjects
    i_male = find(cov.sex==1);i_female = find(cov.sex==0);
    ind_male = randperm(length(i_male));
    ind_female = randperm(length(i_female));
    ind_used = [i_male(ind_male(1:samplesize_th/2));i_female(ind_female(1:samplesize_th/2))];
    
    FC = FC_UKB(ind_used,:);Gender = [cov.sex(ind_used)];Age = [cov.age(ind_used)];
    
    % regressing out age effects for UK Biobank samples
    n_ROI = 94;
    n_person = length(ind_used);
    regresscov = [Age,Age.^2,Age.^3];
    fprintf('starting regressing out covariates\n')

    for i = 1:(n_ROI*(n_ROI-1)/2)
        glmstruct = fitglm(regresscov,FC(:,i));
        b = glmstruct.Coefficients.Estimate;
        FC_regressed(:,i) = FC(:,i)-[ones(n_person,1),regresscov]*b;
    end
    
    fprintf('finish regressing out covariants (UKB)\n')
    
    % preparation for 10-fold CV
    n_person = size(FC_regressed,1);nfold = 10;
    Indices = crossvalind('Kfold',n_person,nfold);
    Output = nan(n_person,1);
    
    for turns = 1:nfold    
        
        % setting training and validation sets according to the assigned label
        Xtrain = FC_regressed(Indices~=turns,:);        Ytrain = Gender(Indices~=turns);
        Xtest = FC_regressed(Indices==turns,:);

        % fitting SVM model. kernelscale=12 is an empirical parameter
        svmstruct = fitcsvm(Xtrain,Ytrain,'Standardize',false,'KernelFunction','linear','KernelScale',12);    
        Output(Indices==turns) = predict(svmstruct,Xtest);    
        fprintf('finish cv %.0f for sample size = %.0f\n',turns,samplesize_th);
        
    end
    
    % calculating 10-fold CV accuracy
    acc(:,rpt) = 1-sum(abs(Gender-Output))/n_person;
    acc_m(:,rpt) = 1-sum(abs(Gender(Gender==1)-Output(Gender==1)))/sum(Gender==1);
    acc_f(:,rpt) = 1-sum(abs(Gender(Gender==0)-Output(Gender==0)))/sum(Gender==0);

    % training on UK Biobank and testing on HCP

    svmstruct = fitcsvm(FC_regressed,Gender,'Standardize',false,'KernelFunction','linear','KernelScale',12);

    fprintf('finish fitting\n')

    % predicting on 4 HCP runs and average the probability to be the final output

    [~,prob1(:,:)] = predict(svmstruct,FC1_regressed);    [~,prob2(:,:)] = predict(svmstruct,FC2_regressed);
    [~,prob3(:,:)] = predict(svmstruct,FC3_regressed);    [~,prob4(:,:)] = predict(svmstruct,FC4_regressed);
    prb1 = normcdf(prob1(:,2));    prb2 = normcdf(prob2(:,2));    prb3 = normcdf(prob3(:,2));    prb4 = normcdf(prob4(:,2));
    gc = mean([prb1,prb2,prb3,prb4],2);
    output_hcp = nan(719,1);    output_hcp(gc>0.5) = 1;    output_hcp(gc<=0.5) = 0;

    % calculating model test accuracy
    acc_hcp(:,rpt) = 1-sum(abs(Gender_hcp-output_hcp))/length(Gender_hcp);
    acc_m_hcp(:,rpt) = 1-sum(abs(Gender_hcp(Gender_hcp==1)-output_hcp(Gender_hcp==1)))/sum(Gender_hcp==1);
    acc_f_hcp(:,rpt) = 1-sum(abs(Gender_hcp(Gender_hcp==0)-output_hcp(Gender_hcp==0)))/sum(Gender_hcp==0);
    fprintf('finish calculating for sample size = %.0f repeating %.0f\n',samplesize_th,rpt);
end