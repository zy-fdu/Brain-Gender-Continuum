% creating permutation set
addpath('permutation functions')
cd data
EB = hcp2blocks('perm_fam.csv');                        % creating exchangeable blocks
[Pset,~] = palm_quickperms(EB,EB,100000);               % creating permutation sets (permute for 100,000 times)
save('Pset.mat','Pset')
cd ..
%%
clc;clear;clf;close;
load('data\gender_spectrum')
% correlation analysis
for i_perm = 1:size(Pset,2)
    % covariates: 2:age, 3:sex, 5:mean FD, 7:handedness
    % calculating linear coefficients
    r_linear(i_perm,:) = partialcorr(gc,pheno{Pset(:,i_perm),:},cov_HCP{:,[2,3,5,7]},'rows','pairwise');
    rm_linear(i_perm,:) = partialcorr(gc(cov_HCP.sex==1),pheno{Pset(cov_HCP.sex==1,i_perm),:},cov_HCP{cov_HCP.sex==1,[2,5,7]},'rows','pairwise');
    rf_linear(i_perm,:) = partialcorr(gc(cov_HCP.sex==0),pheno{Pset(cov_HCP.sex==0,i_perm),:},cov_HCP{cov_HCP.sex==0,[2,5,7]},'rows','pairwise');

    % calculating quadratic coefficients
    r_quad(i_perm,:) = partialcorr(gc.^2,pheno{Pset(:,i_perm),:},[gc,cov_HCP{:,[2,3,5,7]}],'rows','pairwise');
    rm_quad(i_perm,:) = partialcorr(gc(cov_HCP.sex==1).^2,pheno{Pset(cov_HCP.sex==1,i_perm),:},[gc(cov_HCP.sex==1),cov_HCP{cov_HCP.sex==1,[2,5,7]}],'rows','pairwise');
    rf_quad(i_perm,:) = partialcorr(gc(cov_HCP.sex==0).^2,pheno{Pset(cov_HCP.sex==0,i_perm),:},[gc(cov_HCP.sex==0),cov_HCP{cov_HCP.sex==0,[2,5,7]}],'rows','pairwise');

    fprintf('permutation round %.0f\n',i_perm);
end

% calculating p-value
pm = sum(abs(rm)>abs(rm(1,:)))/100000;
pf = sum(abs(rf)>abs(rf(1,:)))/100000;
pm_quad = sum(abs(rm_quad)>abs(rm_quad(1,:)))/100000;
pf_quad = sum(abs(rf_quad)>abs(rf_quad(1,:)))/100000;