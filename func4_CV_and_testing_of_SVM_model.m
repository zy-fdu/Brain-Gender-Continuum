% 10-fold Cross-Validation on UK Biobank
clear;clc;clf;close;
load('data\AAL2_regressedFC.mat')   % the regressed FC was obtained from func3_model_training_preparation function.

% giving labels to each sample, thus setting training and validation sets
nfold = 10;
n_person_UKB = size(FC_regressed,1);n_person_HCP = size(FC1_regressed,1);
Indices = crossvalind('Kfold',n_person_UKB,nfold);
Output = nan(n_person_HCP,1);   % initializing the output vector

% 10-fold cross validation
for turns = 1:nfold
    % setting training and validation sets according to the assigned label
    Xtrain = FC_regressed(Indices~=turns,:);    Ytrain = cov_UKB.sex(Indices~=turns);
    Xtest = FC_regressed(Indices==turns,:);
    
    % fitting SVM model. kernelscale=12 is an empirical parameter
    svmstruct = fitcsvm(Xtrain,Ytrain,'Standardize',false,'KernelFunction','linear','KernelScale',12);
    
    % predicted label and probability(Z-value) were given by the SVM
    [Output(Indices==turns),prob(Indices==turns,:)] = predict(svmstruct,Xtest);
    fprintf('finish cv %.0f\n',turns);

end

% calculating CV accuracy
acc_UKB = 1-sum(abs(cov_UKB.sex-Output))/n_person_UKB;
gc_UKB = normcdf(prob(:,2));    % gender continuum of UKB

% training on UK Biobank and testing on HCP

svmstruct = fitcsvm(FC_regressed,cov_UKB.sex,'Standardize',false,'KernelFunction','linear','KernelScale',12);
fprintf('finish fitting\n')

% predicting on HCP set (4 HCP runs)

[~,prob1] = predict(svmstruct,FC1_regressed);
[~,prob2] = predict(svmstruct,FC2_regressed);
[~,prob3] = predict(svmstruct,FC3_regressed);
[~,prob4] = predict(svmstruct,FC4_regressed);    
fprintf('finish calculating\n');

% turning Z-value to probability
gc1 = normcdf(prob1(:,2));gc2 = normcdf(prob2(:,2));gc3 = normcdf(prob3(:,2));gc4 = normcdf(prob4(:,2));
gc = mean([gc1,gc2,gc3,gc4],2); % gender continuum of HCP samples

% calculating test accuracy
acc = 1-sum(abs(double(gc>0.5)-cov_HCP.sex))/n_person_HCP