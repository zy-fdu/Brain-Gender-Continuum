% model training preparation
clear;clf;close;clc
load('data\AAL2_FC.mat')    % Functional connectivity based on AAL2 parcellation. Including UKB, HCP, IMAGEN datasets, and their corresponding covariates
                            % cov_XXX file includes ID, age, sex, TIV, meanFD, SNR, and handedness
                            % FC file is a n_subject * 4371 (=93*94/2) matrix

% random selecting subjects in order to keep the balance between female and male subjects
i_female = find(cov_UKB.sex==0);i_male = find(cov_UKB.sex==1);
i_f_rand = randperm(length(i_female));i_m_rand = randperm(length(i_male));
i_female_selected = i_female(i_f_rand(1:3700));i_male_selected = i_male(i_m_rand(1:3700));

sub_ind_UKB = [i_female_selected;i_male_selected];
FC_UKB = FC_UKB(sub_ind_UKB,:);
cov_UKB = cov_UKB(sub_ind_UKB,:);

% age, age.^2, and age.^3 were regressed out from FC before training the model.
n_person_HCP = size(FC1_HCP,1);n_person_UKB = size(FC_UKB,1);
regresscov_HCP = [cov_HCP.age,cov_HCP.age.^2,cov_HCP.age.^3];
regresscov_UKB = [cov_UKB.age,cov_UKB.age.^2,cov_UKB.age.^3];
n_ROI = 94;

fprintf('starting regressing out covariates\n')

for i = 1:(n_ROI*(n_ROI-1)/2)
    % regressing out age effect from 4 HCP runs
    glmstruct = fitglm(regresscov_HCP,FC1_HCP(:,i));
    b = glmstruct.Coefficients.Estimate;
    FC1_regressed(:,i) = FC1_HCP(:,i)-[ones(n_person_HCP,1),regresscov_HCP]*b;
    
    glmstruct = fitglm(regresscov_HCP,FC2_HCP(:,i));
    b = glmstruct.Coefficients.Estimate;
    FC2_regressed(:,i) = FC2_HCP(:,i)-[ones(n_person_HCP,1),regresscov_HCP]*b;
    
    glmstruct = fitglm(regresscov_HCP,FC3_HCP(:,i));
    b = glmstruct.Coefficients.Estimate;
    FC3_regressed(:,i) = FC3_HCP(:,i)-[ones(n_person_HCP,1),regresscov_HCP]*b;
    
    glmstruct = fitglm(regresscov_HCP,FC4_HCP(:,i));
    b = glmstruct.Coefficients.Estimate;
    FC4_regressed(:,i) = FC4_HCP(:,i)-[ones(n_person_HCP,1),regresscov_HCP]*b;
    
    % regressing out age effect from UKB
    glmstruct = fitglm(regresscov_UKB,FC_UKB(:,i));
    b = glmstruct.Coefficients.Estimate;
    FC_regressed(:,i) = FC_UKB(:,i)-[ones(n_person_UKB,1),regresscov_UKB]*b;
end


fprintf('finish regressing out covariates.\n')

