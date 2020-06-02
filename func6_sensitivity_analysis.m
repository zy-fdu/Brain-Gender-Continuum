clear;clf;close;clc;
fprintf('loading data...\n')
load('data\AAL2_FC.mat','FC_UKB');
load('data\AAL2_regressedFC.mat')% the regressed FC was obtained from func3_model_training_preparation function. regressed FC of HCP samples were used in this function

addpath('sensitivity analysis functions') % functions for sensitivity analysis were stored in this folder

% preparation for sensitivity analysis
samplesize_th = [250;500;750;1000;1250;1500;1750;2000;2500;3000;3500;4000;5000;6000;7000]; % manually chosen threshold of sample size
delete(gcp('nocreate')) % initializing parallel computing tools
partool = parpool('local',20);
partool.IdleTimeout = Inf;

% model training
parfor i_ssth = 1:length(samplesize_th)
    [acc(i_ssth,:),acc_m(i_ssth,:),acc_f(i_ssth,:),acc_hcp(i_ssth,:),acc_m_hcp(i_ssth,:),acc_f_hcp(i_ssth,:)] = calc_ss(samplesize_th(i_ssth),FC_UKB,cov_UKB,FC1_regressed,FC2_regressed,FC3_regressed,FC4_regressed,cov_HCP.sex);
end

save('results_ss.mat','acc*','samplesize_th')

%% changing age configuration of training set and test its influence to SVM performance
clear;clf;close;clc;
fprintf('loading data...\n')
load('data\AAL2_FC.mat','FC_UKB');
load('data\AAL2_regressedFC.mat') % the regressed FC was obtained from func3_model_training_preparation function. regressed FC of HCP samples were used in this function

addpath('sensitivity analysis functions') % functions for sensitivity analysis were stored in this folder

% processing UKB
choosed_age = [0:300:3000];  % training sample size=3,000 were set, 0%-100% subject older than 65 years old were included. 10% as the step length.

delete(gcp('nocreate'))  % initializing parallel computing tools
partool = parpool('local',20);
partool.IdleTimeout = Inf;

parfor i_percentage_65 = 1:length(choosed_age)
    [acc(i_percentage_65,:),acc_m(i_percentage_65,:),acc_f(i_percentage_65,:),acc_hcp(i_percentage_65,:),acc_m_hcp(i_percentage_65,:),acc_f_hcp(i_percentage_65,:)] = calc65_unreg(choosed_age,i_percentage_65,FC_UKB,cov_UKB,FC1_regressed,FC2_regressed,FC3_regressed,FC4_regressed,cov_HCP.sex);
end

save('results_age_configuration.mat','acc*')

%% changing age upper bound of training samples and test its influence to SVM performance
clear;clf;close;clc;
fprintf('loading data...\n')
load('data\AAL2_FC.mat','FC_UKB','FC1_HCP','FC2_HCP','FC3_HCP','FC4_HCP');  % unregressed FC, used to study age effect
load('data\AAL2_regressedFC.mat'); % the regressed FC was obtained from func3_model_training_preparation function. regressed FC of HCP samples were used in this function

addpath('sensitivity analysis functions') % functions for sensitivity analysis were stored in this folder

age_upperbound = 55:75;
delete(gcp('nocreate'))  % initializing parallel computing tools
partool = parpool('local',20);
partool.IdleTimeout = Inf;

% age and its higher order terms were regressed out in this model
parfor i_ageub = 1:length(age_upperbound)
    [acc(i_ageub,:),acc_m(i_ageub,:),acc_f(i_ageub,:),acc_hcp(i_ageub,:),acc_m_hcp(i_ageub,:),acc_f_hcp(i_ageub,:),ss(i_ageub)] = calc_age_reg(age_upperbound(i_ageub),FC_UKB,cov_UKB,FC1_regressed,FC2_regressed,FC3_regressed,FC4_regressed,cov_HCP.sex);
end
save('results_age_reg.mat','acc*','ss','age_upperbound')
clear acc*

% age and its higher order terms were NOT regressed out in this model
parfor i_ageub = 1:length(age_upperbound)
    [acc(i_ageub,:),acc_m(i_ageub,:),acc_f(i_ageub,:),acc_hcp(i_ageub,:),acc_m_hcp(i_ageub,:),acc_f_hcp(i_ageub,:)] = calc_age_unreg(age_upperbound(i_ageub),FC_UKB,cov_UKB,FC1_HCP,FC2_HCP,FC3_HCP,FC4_HCP,cov_HCP.sex);
end
save('results_unreg_age.mat','acc*','ss','age_upperbound')

%% Excluding brain networks in the model, and test its influence to the performance of SVM
clear;clf;close;clc;
fprintf('loading data...\n')
load('data\AAL2_TC.mat'); % using time-course data to facilitate excluding brain regions from a certain brain network
load('AAL2_network.mat') % file includes 11 brain network and the index of their corresponding brain regions

n_ROI = 94;

% random selecting subjects in order to keep the balance between female and male subjects
i_female = find(cov_UKB.sex==0);i_male = find(cov_UKB.sex==1);
i_f_rand = randperm(length(i_female));i_m_rand = randperm(length(i_male));
i_female_selected = i_female(i_f_rand(1:3700));i_male_selected = i_male(i_m_rand(1:3700));
sub_ind_UKB = [i_female_selected;i_male_selected];

TS_UKB = aal2TC_UKB(sub_ind_UKB);cov_UKB = cov_UKB(sub_ind_UKB,:);
n_person_HCP = length(aal2TC_HCP1);n_person_UKB = length(TS_UKB);

for i_excluded_net = 1:length(index_network)
    clear FC1;clear FC2;clear FC3;clear FC4;clear FC;
    clear FC1_regressed;clear FC2_regressed;clear FC3_regressed;clear FC4_regressed;clear FC_regressed;

    fprintf(['Excluding network %.0f :' network_name{i_excluded_net} '\n'],i_excluded_net)
    % exclude one brain network at a time
    calc_roi = setdiff([1:n_ROI],index_network{i_excluded_net});
    fprintf('remaining %.0f ROIs\n',length(calc_roi))
    
    % calculating FC on 4 HCP runs
    for person = 1:n_person_HCP
        X1 = corr(aal_ts1{person}(:,calc_roi));    X2 = corr(aal_ts2{person}(:,calc_roi));
        X3 = corr(aal_ts3{person}(:,calc_roi));    X4 = corr(aal_ts4{person}(:,calc_roi));
        xtemp1 = [];xtemp2 = [];xtemp3 = [];xtemp4 = [];
        for node = 1:length(calc_roi)
            xtemp1 = [xtemp1,X1(node,node+1:end)];        xtemp2 = [xtemp2,X2(node,node+1:end)];
            xtemp3 = [xtemp3,X3(node,node+1:end)];        xtemp4 = [xtemp4,X4(node,node+1:end)];
        end
        FC1(person,:) = xtemp1;    FC2(person,:) = xtemp2;    FC3(person,:) = xtemp3;    FC4(person,:) = xtemp4;
    end
    % fisher-Z transforming
    FC1 = 0.5*log((1+FC1)./(1-FC1));FC2 = 0.5*log((1+FC2)./(1-FC2));FC3 = 0.5*log((1+FC3)./(1-FC3));FC4 = 0.5*log((1+FC4)./(1-FC4));

    % regressing out age effects from HCP samples
    regresscov = [cov_HCP.age,cov_HCP.age.^2,cov_HCP.age.^3];
    for i = 1:size(FC1,2)
        
        glmstruct1 = fitglm(regresscov,FC1(:,i));b1 = glmstruct1.Coefficients.Estimate;
        FC1_regressed(:,i) = FC1(:,i)-[ones(n_person_HCP,1),regresscov]*b1;
        
        glmstruct2 = fitglm(regresscov,FC2(:,i));b2 = glmstruct2.Coefficients.Estimate;
        FC2_regressed(:,i) = FC2(:,i)-[ones(n_person_HCP,1),regresscov]*b2;
        
        glmstruct3 = fitglm(regresscov,FC3(:,i));b3 = glmstruct3.Coefficients.Estimate;
        FC3_regressed(:,i) = FC3(:,i)-[ones(n_person_HCP,1),regresscov]*b3;
    
        glmstruct4 = fitglm(regresscov,FC4(:,i));b4 = glmstruct4.Coefficients.Estimate;
        FC4_regressed(:,i) = FC4(:,i)-[ones(n_person_HCP,1),regresscov]*b4;
    end
    fprintf('finish regressing out covariates (HCP)\n')

    % processing UKB (same as HCP procedure)
    for person = 1:n_person_UKB
        X = corr(TS_UKB{person}(:,calc_roi));
        xtemp = [];
        for node = 1:length(calc_roi)
            xtemp = [xtemp,X(node,node+1:end)];
        end
        FC(person,:) = xtemp;
    end
    FC = 0.5*log((1+FC)./(1-FC));
    regresscov = [cov_UKB.age,cov_UKB.age.^2,cov_UKB.age.^3];
    for i = 1:size(FC,2)
        glmstruct = fitglm(regresscov,FC(:,i));b = glmstruct.Coefficients.Estimate;
        FC_regressed(:,i) = FC(:,i)-[ones(n_person_UKB,1),regresscov]*b;
    end
    
    fprintf('finish regressing out covariates (UKB)\n')
    
    % 10-fold Cross Validation assessment
    % giving labels to each sample, thus setting training and validation sets
    nfold = 10;
    Indices = crossvalind('Kfold',n_person_UKB,nfold);
    Output = nan(n_person_UKB,1);
    
    for turns = 1:nfold
        % fitting SVM model. kernelscale=12 is an empirical parameter
        svmstruct = fitcsvm(FC_regressed(Indices~=turns,:),cov_UKB.sex(Indices~=turns),'Standardize',false,'KernelFunction','linear','KernelScale',12);
        Output(Indices==turns) = predict(svmstruct,FC_regressed(Indices==turns,:));  
        fprintf('finish cv %.0f with exclusion of cluster %.0f\n',turns,i_excluded_net);
    end
    
    % calculating CV accuracy
    acc(i_excluded_net) = 1-sum(abs(cov_UKB.sex-Output(:)))/n_person_UKB;
    acc_male(i_excluded_net) = 1-sum(abs(cov_UKB.sex(cov_UKB.sex==1)-Output(cov_UKB.sex==1)))/sum(cov_UKB.sex==1);
    acc_female(i_excluded_net) = 1-sum(abs(cov_UKB.sex(cov_UKB.sex==0)-Output(cov_UKB.sex==0)))/sum(cov_UKB.sex==0);
    
    % fitting SVM model and test on HCP dataset
    svmstruct = fitcsvm(FC_regressed,cov_UKB.sex,'Standardize',false,'KernelFunction','linear','KernelScale',12);
    fprintf('finish fitting\n')
    
    % predicting on HCP set
    [~,prob1] = predict(svmstruct,FC1_regressed);[~,prob2] = predict(svmstruct,FC2_regressed);
    [~,prob3] = predict(svmstruct,FC3_regressed);[~,prob4] = predict(svmstruct,FC4_regressed);
    
    % turning Z-value to probability
    gc1 = normcdf(prob1(:,2));gs2 = normcdf(prob2(:,2));gs3 = normcdf(prob3(:,2));gs4 = normcdf(prob4(:,2));
    gc(:,i_excluded_net) = mean([gs1,gs2,gs3,gs4],2);
    fprintf('finish calculating excluding cluster %.0f\n',i_excluded_net);
    
    % calculating test accuracy
    acc_hcp(i_excluded_net) = 1-sum(abs(double(gc(:,i_excluded_net)>0.5)-cov_HCP.sex))/n_person_HCP;
    acc_male_hcp(i_excluded_net) = 1-sum(abs(cov_HCP.sex(cov_HCP.sex==1)-double(gc(cov_HCP.sex==1,i_excluded_net)>0.5)))/sum(cov_HCP.sex==1);
    acc_female_hcp(i_excluded_net) = 1-sum(abs(cov_HCP.sex(cov_HCP.sex==0)-double(gc(cov_HCP.sex==0,i_excluded_net)>0.5)))/sum(cov_HCP.sex==0);

end
% bootstraping for accuracy decline
gc_exc = gc;
load('HCP_SVMresults.mat','gc') % gc obtained func4_CV_and_testing_of_SVM_model function
output_exc = double(gc_exc>0.5);output = double(gc>0.5); % turning gc into output labels
% separate female and male samples to different variables to facilitate further computation
output_f = output(sex_HCP.sex==1);output_exc_f = output_exc(sex_HCP.sex==1,:);
output_m = output(sex_HCP.sex==1);output_exc_m = output_exc(sex_HCP.sex==1,:);
n_f = length(opt_f);n_m = length(opt_m);
for net = 1:11
    for i_bstp = 1:100000
        ind_bstp = randsample(n_person_HCP,n_person_HCP,1); % randomsample of all subjects, female subjects, male subjects
        ind_f = randsample(n_f,n_f,1);ind_m = randsample(n_m,n_m,1);
        
        acc_fullmodel = 1-sum(abs(output(ind_bstp)-cov_HCP.sex(ind_bstp)))/n_person_HCP; % boostrapping in all subjects
        acc_excmodel = 1-sum(abs(output_exc(ind_bstp)-cov_HCP.sex(ind_bstp)))/n_person_HCP;
        acc_decline(i_bstp,net) = acc_fullmodel-acc_excmodel;
        
        acc_f_fullmodel = 1-sum(abs(output_f(ind_f)-0))/n_f;    % boostrapping in females
        acc_f_excmodel = 1-sum(abs(output_exc(ind_f)-0))/n_f;
        acc_f_decline(i_bstp,net) = acc_f_fullmodel-acc_f_excmodel;
        
        acc_m_fullmodel = 1-sum(abs(output(ind_m)-1))/n_m;      % boostrapping in males
        acc_m_excmodel = 1-sum(abs(output_exc(ind_m)-1))/n_m;
        acc_m_decline(i_bstp,net) = acc_m_fullmodel-acc_m_excmodel;
    end
end