%% local_metacognitive_efficiency_analysis
% This script runs the hierchical Bayesian estimates of group metacognitive ratings for local
% confidence, and regressions using the Hierachical model
% in the manuscript, these are presented in Table 2 and Figure 2

% It generates the outputs found on github, which are then used in plotting in the other scripts.

% The hierarachical meta-d' toolbox is located at https://github.com/metacoglab/HMeta-d
% Further information about using the hierarchical model is found at 
% https://github.com/smfleming/HMM/wiki/HMeta-d-tutorial
% It requires JAGS to be on the pathway

load ParticDemogs_and_globals % get demographics data
partics=Partics_and_globals; % and rename it
age_single = partics.age_single';
age_group = partics.age_group';

load trials_mem
load trials_perc

%% First, group metacognitive efficiencies 
for domain =1:2
% Do memory first, then perception
% First, get the organised confidence ratings generated previously by using
% trials2counts during preprocessing_and_exclusions 
    clear nR_S1
    clear nR_S2
 if domain ==1
    nR_S1 = nR_S1_mem; % Choose memory data 
    nR_S2 = nR_S2_mem;
 else 
    nR_S1 = nR_S1_perc; % Choose perception data 
    nR_S2 = nR_S2_perc;
 end
    
%% Fit group data all at once
    clear fit
mcmc_params = fit_meta_d_params; % Choose the default parameters
fit = fit_meta_d_mcmc_group(nR_S1, nR_S2, mcmc_params);
total_HDI=calc_HDI(exp(fit.mcmc.samples.mu_logMratio(:)));
if domain ==1
    fit_mem_all = fit;
    % save ('mem_fit_allsubjects', 'fit_mem_all'); % save if you want
    % save ('mem_allsubjects_HDI', 'total_HDI')
else
    fit_perc_all = fit;
    % save ('perc_fit_allsubjects', 'fit_perc_all');
    % save ('perc_allsubjects_HDI', 'total_HDI')
end

%% Now generate group metacognitive efficiencies for each of the 6 age groups taken singly 

clear eff_6groups_fit 
clear eff_6groups_HDI 
clear eff_6groups_mu 

eff_6groups_fit = []; % create a  variable to put them in
eff_6groups_HDI = [];
eff_6groups_mu = [];

kk = 1;
      nR_S1_singleage = nR_S1(age_group == kk); %find the packagaed confidence ratings
      nR_S2_singleage = nR_S2(age_group == kk);
% Fit model for that group
fit = fit_meta_d_mcmc_group(nR_S1_singleage, nR_S2_singleage, mcmc_params);
eff_6groups_fit.fit1 = fit;
 eff_6groups_mu(kk)= exp(fit.mu_logMratio);
 eff_6groups_HDI(kk,:)=calc_HDI(exp(fit.mcmc.samples.mu_logMratio(:)));
clear fit

kk = 2;
      nR_S1_singleage = nR_S1(age_group == kk); %find the packagaed confidence ratings
      nR_S2_singleage = nR_S2(age_group == kk);
% Fit model for that group
fit = fit_meta_d_mcmc_group(nR_S1_singleage, nR_S2_singleage, mcmc_params);
eff_6groups_fit.fit2 = fit;
 eff_6groups_mu(kk)= exp(fit.mu_logMratio);
 eff_6groups_HDI(kk,:)=calc_HDI(exp(fit.mcmc.samples.mu_logMratio(:)));
clear fit

kk = 3;
      nR_S1_singleage = nR_S1(age_group == kk); %find the packagaed confidence ratings
      nR_S2_singleage = nR_S2(age_group == kk);
% Fit model for that group
fit = fit_meta_d_mcmc_group(nR_S1_singleage, nR_S2_singleage, mcmc_params);
eff_6groups_fit.fit3 = fit;
 eff_6groups_mu(kk)= exp(fit.mu_logMratio);
 eff_6groups_HDI(kk,:)=calc_HDI(exp(fit.mcmc.samples.mu_logMratio(:)));
clear fit

kk = 4;
      nR_S1_singleage = nR_S1(age_group == kk); %find the packagaed confidence ratings
      nR_S2_singleage = nR_S2(age_group == kk);
% Fit model for that group
fit = fit_meta_d_mcmc_group(nR_S1_singleage, nR_S2_singleage, mcmc_params);
eff_6groups_fit.fit4 = fit;
 eff_6groups_mu(kk)= exp(fit.mu_logMratio);
 eff_6groups_HDI(kk,:)=calc_HDI(exp(fit.mcmc.samples.mu_logMratio(:)));
clear fit 

kk = 5;
      nR_S1_singleage = nR_S1(age_group == kk); %find the packagaed confidence ratings
      nR_S2_singleage = nR_S2(age_group == kk);
% Fit model for that group
fit = fit_meta_d_mcmc_group(nR_S1_singleage, nR_S2_singleage, mcmc_params);
eff_6groups_fit.fit5 = fit;
 eff_6groups_mu(kk)= exp(fit.mu_logMratio);
 eff_6groups_HDI(kk,:)=calc_HDI(exp(fit.mcmc.samples.mu_logMratio(:)));
clear fit

kk = 6;
      nR_S1_singleage = nR_S1(age_group == kk); %find the packagaed confidence ratings
      nR_S2_singleage = nR_S2(age_group == kk);
% Fit model for that group
fit = fit_meta_d_mcmc_group(nR_S1_singleage, nR_S2_singleage, mcmc_params);
eff_6groups_fit.fit6 = fit;
 eff_6groups_mu(kk)= exp(fit.mu_logMratio);
 eff_6groups_HDI(kk,:)=calc_HDI(exp(fit.mcmc.samples.mu_logMratio(:)));
clear fit

% Now package outputs
if domain ==1
    mem_6groups_eff.fits = eff_6groups_fit;
    mem_6groups_eff.mu_x6 = eff_6groups_mu;
    mem_6groups_eff.HDI_x6 = eff_6groups_HDI;
else
    perc_6groups_eff.fits = eff_6groups_fit;
    perc_6groups_eff.mu_x6 = eff_6groups_mu;
    perc_6groups_eff.HDI_x6 = eff_6groups_HDI;
end

clear nR_S1_singleage
clear nR_S2_singleage
clear fit
   
end
clear domain
%% Next perform regressions wihtin the extended version of hierachical model of local metacognitive
% efficiency, using fit_meta_d_mcmc_regression
% as described in Hamilton et al (2021)

load('outputs_memory.mat') % Get the possible predictors from memory and perception performance
load('outputs_perception.mat')

for domain = 1:2 % First memory then perception
    clear nR_S1
    clear nR_S2
if domain == 1 % Choose memory data  
   nR_S1 = nR_S1_mem; 
   nR_S2 = nR_S2_mem;
% Prepare covariance, with age, quadratic in age, mean level on the
% diffiulty staircase, and standard deviation of dificulty staircase level
% within each individual
clear cov;
cov = vertcat (normalize(age_single),...
    normalize(age_single.^2),...
    normalize(mem_output_variables.difflevel),...
    normalize(mem_output_variables.diffstd));

else
    nR_S1 = nR_S1_perc; 
    nR_S2 = nR_S2_perc;
clear cov;
cov = vertcat (normalize(age_single),...
    normalize(age_single.^2),...
    normalize(perc_output_variables.difflevel),...
    normalize(perc_output_variables.diffstd));
end

% Fit regression to whole-group data
mcmc_params = fit_meta_d_params; % Choose the default parameters
mcmc_params.estimate_dprime = 0;
fit = fit_meta_d_mcmc_regression(nR_S1, nR_S2, cov, mcmc_params);

HRegression.fit = fit;
HRegression.beta_hdi(1,:) = calc_HDI(fit.mcmc.samples.mu_beta1(:));
HRegression.beta_hdi(2,:) = calc_HDI(fit.mcmc.samples.mu_beta2(:));
HRegression.beta_hdi(3,:) = calc_HDI(fit.mcmc.samples.mu_beta3(:));
HRegression.beta_hdi(4,:) = calc_HDI(fit.mcmc.samples.mu_beta4(:));

if domain ==1
    mem_HRegression = HRegression;
    % file_name = 'mem_HRegression_efficiency'; % for saving if wanted
    % save (file_name,'mem_HRegression' )
else
    perc_HRegression = HRegression;
    % file_name = 'perc_HRegression_efficiency';
    % save (file_name,'perc_HRegression' )
end   
clear file_name
clear nR_S1
clear nR_S2
clear fit
clear HRegression
clear cov
end
clear domain

%% Next, estimate the correlation covariance between memory and 
% perception within individuals using the hierarchical model

% First put into the format for use with fit_meta_d_mcmc_groupCorr

        nR_S1(1).counts = nR_S1_mem ;  
        nR_S2(1).counts = nR_S2_mem ;
       
        nR_S1(2).counts = nR_S1_perc;
        nR_S2(2).counts = nR_S2_perc;

        fit = fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2, mcmc_params);
        corr_estimates_all =[];
        corr_estimates_all.fit = fit;
        corr_estimates_all.hdi_rho(:) = calc_HDI(fit.mcmc.samples.rho(:));
        
%% Fit group data for each of the 6 age groups
corr_estimates_byage =[];
rho_mean = [];
hdi_rho = [];
mcmc_params = fit_meta_d_params;

kk = 1;
clear nR_S1
clear nR_S2
    nR_S1(1).counts = nR_S1_mem(age_group == kk);
    nR_S2(1).counts = nR_S2_mem(age_group == kk);
    nR_S1(2).counts = nR_S1_perc(age_group == kk);
    nR_S2(2).counts = nR_S2_perc(age_group == kk);
fit = fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2, mcmc_params);
rho_mean(kk) = fit.rho;
hdi_rho(kk,:) = calc_HDI(fit.mcmc.samples.rho(:)); 
corr_estimates_byage.fit1 = fit ;
clear fit

kk = 2;
clear nR_S1
clear nR_S2
    nR_S1(1).counts = nR_S1_mem(age_group == kk);
    nR_S2(1).counts = nR_S2_mem(age_group == kk);
    nR_S1(2).counts = nR_S1_perc(age_group == kk);
    nR_S2(2).counts = nR_S2_perc(age_group == kk);
fit = fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2, mcmc_params);
rho_mean(kk) = fit.rho;
hdi_rho(kk,:) = calc_HDI(fit.mcmc.samples.rho(:)); 
corr_estimates_byage.fit2 = fit ;
clear fit

kk = 3;
clear nR_S1
clear nR_S2
    nR_S1(1).counts = nR_S1_mem(age_group == kk);
    nR_S2(1).counts = nR_S2_mem(age_group == kk);
    nR_S1(2).counts = nR_S1_perc(age_group == kk);
    nR_S2(2).counts = nR_S2_perc(age_group == kk);
fit = fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2, mcmc_params);
rho_mean(kk) = fit.rho;
hdi_rho(kk,:) = calc_HDI(fit.mcmc.samples.rho(:)); 
corr_estimates_byage.fit3 = fit ;
clear fit

kk = 4;
clear nR_S1
clear nR_S2
    nR_S1(1).counts = nR_S1_mem(age_group == kk);
    nR_S2(1).counts = nR_S2_mem(age_group == kk);
    nR_S1(2).counts = nR_S1_perc(age_group == kk);
    nR_S2(2).counts = nR_S2_perc(age_group == kk);
fit = fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2, mcmc_params);
rho_mean(kk) = fit.rho;
hdi_rho(kk,:) = calc_HDI(fit.mcmc.samples.rho(:)); 
corr_estimates_byage.fit4 = fit ;
clear fit

kk = 5;
clear nR_S1
clear nR_S2
    nR_S1(1).counts = nR_S1_mem(age_group == kk);
    nR_S2(1).counts = nR_S2_mem(age_group == kk);
    nR_S1(2).counts = nR_S1_perc(age_group == kk);
    nR_S2(2).counts = nR_S2_perc(age_group == kk);
fit = fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2, mcmc_params);
rho_mean(kk) = fit.rho;
hdi_rho(kk,:) = calc_HDI(fit.mcmc.samples.rho(:)); 
corr_estimates_byage.fit5 = fit ;
clear fit

kk = 6;
clear nR_S1
clear nR_S2
    nR_S1(1).counts = nR_S1_mem(age_group == kk);
    nR_S2(1).counts = nR_S2_mem(age_group == kk);
    nR_S1(2).counts = nR_S1_perc(age_group == kk);
    nR_S2(2).counts = nR_S2_perc(age_group == kk);
fit = fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2, mcmc_params);
rho_mean(kk) = fit.rho;
hdi_rho(kk,:) = calc_HDI(fit.mcmc.samples.rho(:)); 
corr_estimates_byage.fit6 = fit ;
clear fit
clear kk
clear nR_S1
clear nR_S2

corr_estimates_byage.rho = rho_mean;
corr_estimates_byage.HDI = hdi_rho;
