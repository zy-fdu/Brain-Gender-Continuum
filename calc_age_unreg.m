function [acc,acc_m,acc_f,acc_hcp,acc_m_hcp,acc_f_hcp] = calc_age_unreg(age_ul,FC_UKB,cov,FC1,FC2,FC3,FC4,Gender_hcp)
for rpt = 1:100
    
    % random selecting subjects in order to keep the balance between female and male subjects
    i_male = find(cov.sex==1 & cov.age<=age_ul);
    i_female = find(cov.sex==0 & cov.age<=age_ul);
    ind_male = randperm(length(i_male));
    ind_female = randperm(length(i_female));
    ss = min(length(i_male),length(i_female))*2;
    ind_use = [i_male(ind_male(1:ss/2));i_female(ind_female(1:ss/2))];
    
    FC = FC_UKB(ind_use,:);    Gender = [cov.sex(ind_use)];    Age = [cov.age(ind_use)];
    
    % preparation for 10-fold CV
    n_person = size(FC,1);    nfold = 10;    Indices = crossvalind('Kfold',n_person,nfold);
    Output = nan(n_person,1);
    
    for turns = 1:nfold    
        % setting training and validation sets according to the assigned label
        Xtrain = FC(Indices~=turns,:);        Ytrain = Gender(Indices~=turns);
        Xtest = FC(Indices==turns,:);
        
        % fitting SVM model. kernelscale=12 is an empirical parameter
        svmstruct = fitcsvm(Xtrain,Ytrain,'Standardize',false,'KernelFunction','linear','KernelScale',12);    
        Output(Indices==turns) = predict(svmstruct,Xtest);    
        fprintf('finish cv %.0f for upper bound of age = %.0f\n',turns,age_ul);
        
    end
    
    % calculating CV accuracy
    acc(rpt) = 1-sum(abs(Gender-Output))/n_person;
    acc_m(rpt) = 1-sum(abs(Gender(Gender==1)-Output(Gender==1)))/sum(Gender==1);
    acc_f(rpt) = 1-sum(abs(Gender(Gender==0)-Output(Gender==0)))/sum(Gender==0);

    % training on UK Biobank and testing on HCP
    svmstruct = fitcsvm(FC,Gender,'Standardize',false,'KernelFunction','linear','KernelScale',12);
    fprintf('finish fitting\n')

    % predicting on 4 HCP runs and average the probability to be the final output
    [~,prob1(:,:)] = predict(svmstruct,FC1);    [~,prob2(:,:)] = predict(svmstruct,FC2);
    [~,prob3(:,:)] = predict(svmstruct,FC3);    [~,prob4(:,:)] = predict(svmstruct,FC4);
    prb1 = normcdf(prob1(:,2));    prb2 = normcdf(prob2(:,2));    prb3 = normcdf(prob3(:,2));    prb4 = normcdf(prob4(:,2));
    gc = mean([prb1,prb2,prb3,prb4],2);
    output_hcp = nan(719,1);    output_hcp(gc>0.5) = 1;    output_hcp(gc<=0.5) = 0;

    % calculating test accuracy
    acc_hcp(rpt) = 1-sum(abs(Gender_hcp-output_hcp))/length(Gender_hcp);
    acc_m_hcp(rpt) = 1-sum(abs(Gender_hcp(Gender_hcp==1)-output_hcp(Gender_hcp==1)))/sum(Gender_hcp==1);
    acc_f_hcp(rpt) = 1-sum(abs(Gender_hcp(Gender_hcp==0)-output_hcp(Gender_hcp==0)))/sum(Gender_hcp==0);
    fprintf('finish calculating for upper bound of age = %.0f repeating %.0f\n',age_ul,rpt);
end