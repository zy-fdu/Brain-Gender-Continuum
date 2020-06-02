# Brain-continuum
Key words: the UK Biobank cohort; the Human Connectome Project cohort; the IMAGEN cohort (FU2); multi-level block permutation; support vector machine
Data and Code of paper:
XXXXXXXXXXX

NOTE

1.Per the data accessing policy of the GRIP (http://www.gripinfo.ca/Grip/Public/www/QuoiNeuf/en/faitssaillants.asp?langue=en),
the QLSCD data could be used only after approval. Therefore, the QLSCD data in this folder are not allowed to be shared to anyone else.
2.Per the data accessing policy of the NDA agreement, the usage of the ABCD data needs the approval by NIMH. The ABCD cohort Release v1.0 (DOI: 10.15154/1460410) was used in this study.
3.The aal2TC.mat and aal2FC.mat file were time courses and functional connectivity of samples in all three dataset, which were preprocessed by Weikang Gong, who was a colleague in our lab and currently a PhD student at University of Oxford. The pipeline descripted in the Method section were used to preprocess the data.
4.If the data and codes are used in your work, please cite the above reference, namely XXXXX.

SUMMARY （Matlab2019a; PALM:https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM)

step 1: Sex differences in global functional connectivity. (func1_sex_difference_of_gFC.m)

  step 1.1: Fitting localized spline curve for global functional connectivity without considering the covariates.
  
  step 1.2: Fitting localized spline curve considering the covariates.
  
  step 1.3: Comparison of different fitted spline curves.
  
  step 1.4: Validating the findings on gFC using YMU dataset.


step 2: Sex differences in local functional connectivity. (func2_local_sex_difference.m)

  step 2.1: Regional-level comparison of sex differences of functional connectivity.
  
  step 2.2: Network-level comparison of sex differences of functional connectivity.


step 3: Balancing training samples in two sexes and regressing out age and their higher order terms to make preparations for the training. (func3_model_training_preparation.m)

  step 3.1: Randomly choosing equal female and male subjects to keep the balance between two sexes in training samples.
  
  step 3.2: Regressing out age and their higher order terms from functional connectivity from training and test sets respectively.


step 4: Training and assessment of the performance of the SVM model. (func4_CV_and_testing_of_SVM_model.m)

  step 4.1: 10-fold cross validation assessment on UK Biobank dataset.
  
  step 4.2: Test set assessment on HCP dataset.


step 5: Testing the validity of the SVM model to make sure androgyny given by the SVM was between female and male rather than noise. (func5_validity_of_gender_continuum.m)


step 6: Sensitivity analysis of the SVM model. (func6_sensitivity_analysis.m)

  step 6.1: Changing sample size of training set, and test its influence on test accuracy.
  
  step 6.2: Changing age configuration (percentage of training samples older than 65 years old), and test its influence on test accuracy.
  
  step 6.3: Changing the age upper bound of training samples, and test its influence on test accuracy.
  
  
    step 6.3.1: Age and its higher order terms were regressed from the functional connectivity as what has been done in other part of the study.
    
    step 6.3.2: Age and its higher order terms were not regressed from the functional connectivity to study the age effect on test accuracy.
    
  step 6.4: Excluding one brain network used in the model at a time, and test its influence on test accuracy.
  

step 7: Apply gender continuum given by the SVM model to internalizing problems. (func7_application.m)

  step 7.1: define the multi-level blocks according to family structure (perm_fam.csv) in the HCP S900 data, and generate the block permuted sample by calling the PALM function.
  
  step 7.2: correlation analysis (both linearly and quadratically) of the relationship between gender continuum and internalizing problems.
