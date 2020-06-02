% sex differences of gFC
clear;clc;clf;close;
load('data\AAL2_FC.mat')    % Functional connectivity based on AAL2 parcellation. Including UKB, HCP, IMAGEN datasets, and their corresponding covariates
                            % cov_XXX file includes ID, age, sex, TIV, meanFD, SNR, and handedness
                            % FC file is a n_subject * 4371 (=93*94/2) matrix

% concatenation of gFC and age for female subjects
gFC = [globalFC_IMAGEN(cov_IMAGEN.sex==0);globalFC_HCP(cov_HCP.sex==0);globalFC_UKB(cov_UKB.sex==0)];
age = [cov_IMAGEN.age(cov_IMAGEN.sex==0);cov_HCP.age(cov_HCP.sex==0);cov_UKB.age(cov_UKB.sex==0)];

% fit and plot gFC for female subjects
splinemdl_female_nocov = fit(age,gFC,'smoothingspline','SmoothingParam',0.001);
plot(splinemdl_female_nocov,age,gFC)
hold on

% concatenation of gFC and age for male subjects
gFC = [globalFC_IMAGEN(cov_IMAGEN.sex==1);globalFC_HCP(cov_HCP.sex==1);globalFC_UKB(cov_UKB.sex==1)];
age = [cov_IMAGEN.age(cov_IMAGEN.sex==1);cov_HCP.age(cov_HCP.sex==1);cov_UKB.age(cov_UKB.sex==1)];

% fit and plot gFC for male subjects
splinemdl_male_nocov = fit(age,gFC,'smoothingspline','SmoothingParam',0.001);
plot(splinemdl_male_nocov,age,gFC)


%% fit gFC for female subjects with covariates

% concatenation of TIV, meanFD and gFC, sex and age
TIV = [cov_IMAGEN.TIV;cov_HCP.TIV;cov_UKB.TIV];
meanFD = [cov_IMAGEN.meanFD;cov_HCP.meanFD;cov_UKB.meanFD];
gFC = [globalFC_IMAGEN;globalFC_HCP;globalFC_UKB];
sex_all = [cov_IMAGEN.sex;cov_HCP.sex;cov_UKB.sex];
age_all = [cov_IMAGEN.age;cov_HCP.age;cov_UKB.age];

% normalization and concatenation of SNR
SNR = [(cov_IMAGEN.SNR-mean(cov_IMAGEN.SNR))/std(cov_IMAGEN.SNR);(cov_HCP.SNR-mean(cov_HCP.SNR))/std(cov_HCP.SNR);(cov_UKB.SNR-mean(cov_UKB.SNR))/std(cov_UKB.SNR)];

cov_all = [TIV,meanFD,SNR];

% regressing out covariates from gFC
mdl = fitglm(cov_all,gFC);
gFC_res = gFC - [ones(length(sex_all),1),cov_all]*mdl.Coefficients.Estimate;

% fit and plot gFC for female subjects with covariates
splinemdl_female_cov = fit(age_all(sex_all==0),gFC_res(sex_all==0),'smoothingspline','SmoothingParam',0.001);
plot(splinemdl_female_cov,age_all(sex_all==0),gFC_res(sex_all==0))
hold on

% fit and plot gFC for male subjects with covariates
splinemdl_male_cov = fit(age_all(sex_all==1),gFC_res(sex_all==1),'smoothingspline','SmoothingParam',0.001);
plot(splinemdl_male_cov,age_all(sex_all==1),gFC_res(sex_all==1))

%% comparison of different trajectories
x = [17:0.2:78];    % step length of 0.2 years was set in our study

% calculating fitted value for selected age values
y_res_f = splinemdl_female_cov(x);y_res_m = splinemdl_male_cov(x);
y_f = splinemdl_female_nocov(x);y_m = splinemdl_male_nocov(x);

% comparing two trajectories
[r_f,p_f] = corr(y_res_f,y_f);[r_m,p_m] = corr(y_res_m,y_m);

% difference of two trajectories
y_diff_2res = y_res_m - y_res_f;y_diff = y_m - y_f;

%% validation of findings on gFC using YMU dataset
load('data\YMU_data.mat')

plot(YMU_age(YMU_sex==0),globalFC(YMU_sex==0),'r+');
hold on
plot(YMU_age(YMU_sex==1),globalFC(YMU_sex==1),'b+');
lsline
legend('Female','Male')

% boostrap comparing trajectory of two sexes
i_m = find(YMU_sex==1);
i_f = find(YMU_sex==0);
for a = 1:100000    % 100,000 repetitions for bootstrap
    ind_m = randsample(length(i_m),length(i_m),1);
    ind_f = randsample(length(i_f),length(i_f),1);
    r_m(a) = corr(gFC_YMU(i_m(ind_m)),YMU_age(i_m(ind_m)));
    r_f(a) = corr(gFC_YMU(i_f(ind_f)),YMU_age(i_f(ind_f)));
end