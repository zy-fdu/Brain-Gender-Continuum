% local sex differences
clear;clc;clf;close;
load('data\AAL2_FC.mat')    % Functional connectivity based on AAL2 parcellation. Including UKB, HCP, IMAGEN datasets, and their corresponding covariates
                            % cov_XXX file includes ID, age, sex, TIV, meanFD, SNR, and handedness
                            % FC file is a n_subject * 4371 (=93*94/2) matrix

n_ROI = 94;

% edge-wise comparison
for ind_edge = 1:n_ROI*(n_ROI-1)/2
    % d is the vectorized effect size, d_cov is the vectorized effect size controlling covariates
    % fit model for IMAGEN participants (cov 7-13 : 7 scanning sites)
    mdl = fitglm([cov_IMAGEN.sex,cov_IMAGEN{:,7:13}],FC_IMAGEN(:,ind_edge));
    n1 = sum(cov_IMAGEN.sex==0);n2 = sum(cov_IMAGEN.sex==1);
    d(ind_edge,1) = mdl.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2)); % calculate cohen's d as effect size
    
    % fit model for HCP participants
    mdl = fitglm([cov_HCP.sex],FC_HCP(:,ind_edge));
    n1 = sum(cov_HCP.sex==0);n2 = sum(cov_HCP.sex==1);
    d(ind_edge,2) = mdl.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));

    % fit model for UK Biobank participants
    mdl = fitglm([cov_UKB.sex],FC_UKB(:,ind_edge));
    n1 = sum(cov_UKB.sex==0);n2 = sum(cov_UKB.sex==1);
    d(ind_edge,3) = mdl.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));
    
    % fit model for IMAGEN participants controlling covariates
    % (cov 3 : sex; cov 4:TIV; cov 5:mean FD, cov 6:SNR, cov7-13: scanning sites, 3-6 has same meaning for all datasets)
    mdl = fitglm([cov_IMAGEN{:,[3:13]}],FC_IMAGEN(:,ind_edge));
    n1 = sum(cov_IMAGEN.sex==0);n2 = sum(cov_IMAGEN.sex==1);
    d_cov(ind_edge,1) = mdl.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));
    
    % fit model for HCP participants controlling covariates
    mdl = fitglm([cov_HCP{:,[3:6]}],FC_HCP(:,ind_edge));
    n1 = sum(cov_HCP.sex==0);n2 = sum(cov_HCP.sex==1);
    d_cov(ind_edge,2) = mdl.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));
    
    % fit model for UK Biobank participants controlling covariates    
    mdl = fitglm([cov_UKB{:,[3:6]}],FC_UKB(:,ind_edge));
    n1 = sum(cov_UKB.sex==0);n2 = sum(cov_UKB.sex==1);
    d_cov(ind_edge,3) = mdl.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));
    
    fprintf('calculating %.0f FC\n',ind_edge)
    
end

% converting vectorized results into matrices
clear indmat_aal
n = 0;
% constructing index matrix
for i_roi_y = 1:n_ROI
    for i_roi_x = i_roi_y+1:n_ROI
        n = n+1;
        indmat_aal(i_roi_x,i_roi_y) = n;
    end
end
indmat_aal = [indmat_aal,zeros(n_ROI,1)];

D = zeros(n_ROI,n_ROI,3);D_cov = zeros(n_ROI,n_ROI,3);
% D is the effect size matrix, D_cov is the effect size matrix controlling covariates
for i_dataset = 1:3
    % initialization
    temp = zeros(n_ROI,n_ROI);tempcov = zeros(n_ROI,n_ROI);
    % locating position of each edge in the matrix
    for i_edge = 1:size(d,1)
        temp(indmat_aal==i_edge) = d(i_edge,i_dataset);tempcov(indmat_aal==i_edge) = d_cov(i_edge,i_dataset);
    end
    % constructing result matrix
    D(:,:,i_dataset) = temp+temp';D_cov(:,:,i_dataset) = tempcov+tempcov';
end

%
%% network-wise comparison
% net_XXX file is the extracted network-level connectivity. it is a 11*11*n_subject tensor (11 networks)
for net_ind_a = 1:11
    for net_ind_b = 1:11
        % network level analysis on IMAGEN
        glm_nocov = fitglm([cov_IMAGEN.sex,cov_IMAGEN{:,7:13}],squeeze(net_IMAGEN(net_ind_a,net_ind_b,:)));
        n1 = sum(cov_IMAGEN{:,3}==0);n2 = sum(cov_IMAGEN{:,3}==1);
        d_net(net_ind_a,net_ind_b,1) = glm_nocov.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));
        
        % network level analysis on IMAGEN controlling covariates
        glm_cov = fitglm(cov_IMAGEN{:,[3:13]},squeeze(net_IMAGEN(net_ind_a,net_ind_b,:)));
        n1 = sum(cov_IMAGEN{:,3}==0);n2 = sum(cov_IMAGEN{:,3}==1);
        d_cov_net(net_ind_a,net_ind_b,1) = glm_cov.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));
        
        % network level analysis on HCP
        glm_nocov = fitglm(cov_HCP.sex,squeeze(net_HCP(net_ind_a,net_ind_b,:)));
        n1 = sum(cov_HCP{:,3}==0);n2 = sum(cov_HCP{:,3}==1);
        d_net(net_ind_a,net_ind_b,2) = glm_nocov.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));
        
        % network level analysis on HCP controlling covariates
        glm_cov = fitglm(cov_HCP{:,[3:6]},squeeze(net_HCP(net_ind_a,net_ind_b,:)));
        n1 = sum(cov_HCP{:,3}==0);n2 = sum(cov_HCP{:,3}==1);
        d_cov_net(net_ind_a,net_ind_b,2) = glm_cov.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));
        
        % network level analysis on UK Biobank
        glm_nocov = fitglm(cov_UKB.sex,squeeze(net_UKB(net_ind_a,net_ind_b,:)));
        n1 = sum(cov_UKB{:,3}==0);n2 = sum(cov_UKB{:,3}==1);
        d_net(net_ind_a,net_ind_b,3) = glm_nocov.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));
        
        % network level analysis on UK Biobank controlling covariates
        glm_cov = fitglm(cov_UKB{:,[3:6]},squeeze(net_UKB(net_ind_a,net_ind_b,:)));
        n1 = sum(cov_UKB{:,3}==0);n2 = sum(cov_UKB{:,3}==1);
        d_cov_net(net_ind_a,net_ind_b,3) = glm_cov.Coefficients.tStat(2)*(n1+n2)/(sqrt((n1+n2-2)*n1*n2));
    end
end