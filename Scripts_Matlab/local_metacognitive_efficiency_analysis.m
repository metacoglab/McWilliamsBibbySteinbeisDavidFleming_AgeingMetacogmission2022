%% local_metacognitive_efficiency_analysis
% This script performs hierchical Bayesain estimates of group metacognitive ratings for local
% confidence, and regressions using the Hierachical model
% in the manuscript, these are presented in Table 2 and Figure 2

% The hierarachical meta-d' toolbox is located at https://github.com/metacoglab/HMeta-d
% This analysis requires helper functions offit_meta_d_params, fit_meta_d_mcmc_group, calc_HDI 
% and text files Bayes_metad_group.txt
% Further information about using the hierarchical model is found at 
% https://github.com/smfleming/HMM/wiki/HMeta-d-tutorial
% It requires JAGS to be on the pathway

% Needs helper functions of:
% fit_meta_d_mcmc_group
% and text files of:
% fit_meta_d_params

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
    save ('mem_fit_allsubjects', 'fit_mem_all');
    save ('mem_allsubjects_HDI', 'total_HDI')
else
    fit_perc_all = fit;
    save ('perc_fit_allsubjects', 'fit_perc_all');
    save ('perc_allsubjects_HDI', 'total_HDI')
end

%% Generate group metacognitive efficiencies for each of the 6 age groups 

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

% Now re-name and save outputs

if domain ==1
    mem_6groups_eff.fits = eff_6groups_fit;
    mem_6groups_eff.mu_x6 = eff_6groups_mu;
    mem_6groups_eff.HDI_x6 = eff_6groups_HDI;

save ('mem_6group_eff_fits', 'mem_6groups_eff')

else
    perc_6groups_eff.fits = eff_6groups_fit;
    perc_6groups_eff.mu_x6 = eff_6groups_mu;
    perc_6groups_eff.HDI_x6 = eff_6groups_HDI;

save ('perc_6group_eff_fits', 'perc_6groups_eff')
end

clear nR_S1_singleage
clear nR_S2_singleage
clear fit
   
end
clear domain
%% Now perform regressions within the hierachical model of local metacognitive
% efficiency, using fit_meta_d_mcmc_regression
% as described in Hamilton et al (2021)
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
    normalize(memory_variables.difflevel),...
    normalize(memory_variables.diffstd));

else
    nR_S1 = nR_S1_perc; 
    nR_S2 = nR_S2_perc;
clear cov;
cov = vertcat (normalize(age_single),...
    normalize(age_single.^2),...
    normalize(perception_variables.difflevel),...
    normalize(perception_variables.diffstd));
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
    file_name = 'mem_HRegression_efficiency';
    save (file_name,'mem_HRegression' )
else
    perc_HRegression = HRegression;
    file_name = 'perc_HRegression_efficiency';
    save (file_name,'perc_HRegression' )
end   
clear file_name
clear nR_S1
clear nR_S2
clear fit
clear HRegression

end

%% Plotting efficiencies with each age with HDIs
% Figure 2A
% Get colour scheme
col = [[230, 10, 0]/255;... %Red
    [51, 140, 255]/255;... %Blue
    [60,179,113]/255;...%Green
    [255, 204, 204]/255;...%Light Red
    [100, 220, 255] / 255;...%Light Blue
    [152,251,152]/255;... %Light Green
    [85,107,47]/255];%Deep Green

for domain = 1:2
if domain ==1  % choose the cognitive domain
    eff_6groups_HDI = mem_6groups_eff.HDI_x6;
    eff_6groups_mu  = mem_6groups_eff.mu_x6;
    figure(21)

else
    eff_6groups_HDI = perc_6groups_eff.HDI_x6;
    eff_6groups_mu  = perc_6groups_eff.mu_x6;
    figure(22)
end

set(gcf, 'Position', [800 300 350 290]); % 

hold('on');
hold('all');
% Plot HDIs for each of the 6 age groups
for kk = 1:6
line ([age_groupmeans(kk), age_groupmeans(kk)],[eff_6groups_HDI(kk,1), eff_6groups_HDI(kk,2)], 'Color',col(domain+3,:),'LineWidth',3);
end
clear kk 

hold on
% Plot mean estimates
hscat = scatter (age_groupmeans, eff_6groups_mu, 60, 'MarkerEdgeColor', col(domain, :),'MarkerFaceColor', col(domain, :));
hline = line (age_groupmeans, eff_6groups_mu, 'Color', col(domain, :), 'LineWidth',3);
hold on
   
%% Set axes
axh = gca;
axh.FontSize = 14;
axh.FontName = 'Arial';
x_label = 'age group';
x_limits = [17.8 74];
x_ticks = age_groupmeans;
x_ticklabels = [{'18-27'},{'28-37'},{'38-47'},{'48-57'},{'58-67'},{'68+'}];

y_label =  'metacognitive efficiency';
if domain ==1
    y_limits = [0.7, 1.4]; 
    y_ticklabels = 0.7:0.1:1.4;
else
    y_limits = [0.3, 1]; 
    y_ticklabels = 0.3:0.1:1;
end

xlabel(x_label, 'FontSize' , 14)
ylabel(y_label, 'FontSize' , 14)
xlim(x_limits)
ylim(y_limits)
axh.XTick = x_ticks;
axh.XTickLabel = x_ticklabels;
axh.YTick = y_ticklabels;
xtickangle(45)
% axh.YRuler.TickLabelGapOffset = 0; % 4; % move tick labels slightly 

box off
hold('all');

if  domain == 1
    fig_filename = 'Figure2Ai'; 
    savefig (gcf,fig_filename) 
    saveas(gcf,fig_filename, 'pdf') 
    clear fig_filename    

else
    fig_filename = 'Figure2Aii'; 
    savefig (gcf,fig_filename) 
    saveas(gcf,fig_filename, 'pdf') 
    clear fig_filename
end
end

%% Plotting betas with HDIs
% Figure 2B

figure (23)
hold on

mem_beta_hdi = mem_HRegression.beta_hdi;
perc_beta_hdi = perc_HRegression.beta_hdi; 

domain = 1; % first memory
for kk = 1:4    
mem_lines = line ([kk, kk ],[mem_beta_hdi(kk,1), mem_beta_hdi(kk,2)], 'Color',col(domain+3,:),'LineWidth',3);
end
clear kk 
hold on
scatter(1:4, [mem_HRegression.fit.mu_beta1,mem_HRegression.fit.mu_beta2,mem_HRegression.fit.mu_beta3,mem_HRegression.fit.mu_beta4], 60,'MarkerEdgeColor',col(domain,:), 'MarkerFaceColor',col(domain,:));

domain = 2; % then perception
for kk = 1:4
    hold on
 perc_lines = line ([kk+0.3, kk+0.3],[perc_beta_hdi(kk,1), perc_beta_hdi(kk,2)], 'Color',col(domain+3,:),'LineWidth',3);
end
clear kk 
hold on
scatter([1.3, 2.3, 3.3, 4.3], [perc_HRegression.fit.mu_beta1,perc_HRegression.fit.mu_beta2,perc_HRegression.fit.mu_beta3,perc_HRegression.fit.mu_beta4], 60,'MarkerEdgeColor',col(domain,:), 'MarkerFaceColor',col(domain,:));

 x_limits = ([0.8 4.8]);
 y_limits = ([-0.8, 0.85 ]);
x_ticks =  [1.15,2.15, 3.15, 4.15];
 y_ticks = [-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8];
 x_label = 'covariate';
 y_label =  'sampled regression betas';
x_ticklabels = {'age','age squared','difficulty mean','difficulty std'};

set(gca, 'FontSize', 12);
xlabel (x_label)
ylabel (y_label)
xlim (x_limits)
ylim (y_limits)
xticks(x_ticks)
yticks(y_ticks)
xticklabels(x_ticklabels) 

yline (0);
% legend([mem_lines perc_lines],{'memory','perception'},'Location','northeast')

file_name = 'Figure2B_RHMetaD';
savefig (gcf, file_name);
saveas (gcf, file_name, 'pdf');
clear file_name

clear x_limits 
clear y_limits
clear x_ticks 
clear y_ticks 
clear x_label 
clear y_label 
clear x_ticklabels

clear mem_lines
clear perc_lines

%% Next, estimate of the correlation covariance between memory and 
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
        
        save('domaingeneral_efficiency_all', 'corr_estimates_all');
        
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

save ( 'domaingeneral_effic_6ages.mat', 'corr_estimates_byage');

% Plot the 6 age groups dervied from group_Corr
figure(24) 
set(gcf, 'Color',[1 1 1],'Position', [650 100 600 360]) 
% set(gcf, 'Position', [800 400 234 245])
% set(gcf, 'Position', [800 400 350 290]) % for bias and pre-, post-globals
% set(gcf, 'Position', [800 400 420 360]) % for eff

x_label = 'age group';
x_limits = [17.8 75];
x_ticks = age_groupmeans;
x_ticklabels = [{'18-27'},{'28-37'},{'38-47'},{'48-57'},{'58-67'},{'68+'}];
xtickangle(35)

y_limits = [-0.5, 1];
y_ticks = [-0.4 -0.2, 0 0.2 0.4 0.6 0.8 1];
y_ticklabels = [-0.4 -0.2, 0 0.2 0.4 0.6 0.8 1];
yline(0, 'LineWidth', 1);

xlim(x_limits)
ylim(y_limits) 

clear axDG
axDG = gca;
axDG.FontSize = 14;
axDG.FontName = 'Arial';
axDG.XTick = x_ticks;
axDG.YTick = y_ticks;
axDG.XTickLabel = x_ticklabels;
axDG.YTickLabel = y_ticklabels;
axDG.YRuler.TickLabelGapOffset = 0;  % move tick labels slightly left or right

box('off');
hold('all');

%% Scatter plot the data
hold(axDG,'on');

group_HDI = corr_estimates_byage.HDI;
metric_group_means = corr_estimates_byage.rho;

for kk = 1:6
clear line_for_HDIs
line_for_HDIs = line([age_groupmeans(kk),age_groupmeans(kk)],...
    [group_HDI(kk,1),group_HDI(kk,2)]);
    set(line_for_HDIs,'Color',col(3,:),'LineWidth',5)
hold on
end
clear kk

hmet = plot (age_groupmeans, metric_group_means, '-o',...
    'color', col(7,:),...
    'LineWidth',4,...
    'MarkerSize',10,...
    'MarkerEdgeColor',col(7,:),...
    'MarkerFaceColor',col(7,:));
clear line_for_HDIs
clear x_label 
clear x_limits 
clear x_ticks 
clear x_ticklabels 
clear y_limits 
clear y_ticks 
clear y_ticklabels
clear rho_mean

%% Save
fig_filename = 'Fig2C_effDG';
savefig (gcf,fig_filename) 
saveas(gcf,fig_filename, 'pdf') 
clear fig_filename

%% Finally, plot the d-primes against age for Figure 1C, 
% now that the d-primes are available

% Set the parameters and choose the domain
  jj=1;
  while jj <3
x_var = age_single; % plots on age first
x_limits = [17.8 85];
x_ticks = age_groupmeans;
x_ticklabels = [{'18-27'},{'28-37'},{'38-47'},{'48-57'},{'58-67'},{'68+'}];

 if jj == 1
    fig_filename = 'Fig1Ci_mem_dprime_on_age';
    domain = 1; % 1 for mem, 2 for perc
    figure(13) 
    y_var = fit_mem_all.d1; 
    y_limits = [0,3.4];
 elseif jj == 2 
    fig_filename = 'Fig1Cii_perc_dprime_on_age';
    domain = 2; 
    figure(14) 
    y_var = fit_perc_all.d1; 
    y_limits = [0,3.4];
end
 
% Draw figure
set(gcf, 'Position', [800 400 190 245],'Color',[1,1,1]);

box('off');
hold('all');

%% Scatter plot the data
scatter ( x_var, y_var, 12, 'MarkerEdgeColor', col(domain+3, :),'MarkerFaceColor', col(domain+3, :));
hold on
 
%% Add means for the 6 age groups
for kk = 1:6
y_group_means(kk) = mean (y_var(age_group==kk)); % mean of the metric within each age group
y_group_std(kk) = std (y_var(age_group==kk));
end
clear kk

plot (age_groupmeans, y_group_means, '-o',...
    'color', col(domain,:),...
    'LineWidth',4,...
    'MarkerSize',8,...
    'MarkerEdgeColor',col(domain,:),...
    'MarkerFaceColor',col(domain,:));

%% Add s.d.s for the 6 age groups
for kk = 1:6
    clear line_for_std
line_for_std = line([age_groupmeans(kk),age_groupmeans(kk)],...
    [y_group_means(kk)- y_group_std(kk),...
    y_group_means(kk)+ y_group_std(kk)]);
    set(line_for_std,'Color',col(domain,:),'LineWidth',4);
hold on
end
clear kk

%% Set axes
axh = gca;
axh.FontSize = 10;
axh.FontName = 'Arial';
axh.XTick = x_ticks;
axh.XTickLabel = x_ticklabels;
axh.YRuler.TickLabelGapOffset = -1;
% set(axh,'LooseInset',get(axh,'TightInset'));

xtickangle(45)
xlim(x_limits)
ylim(y_limits)

% Tighthen up margins
 tightInset = get(gca, 'TightInset');
position(1) = tightInset(1);
position(2) = tightInset(2);
position(3) = 1 - tightInset(1) - tightInset(3);
position(4) = 1 - tightInset(2) - tightInset(4);
set(axh, 'Position', position);

% now save the figure 
    savefig (gcf,fig_filename) 
    saveas(gcf,fig_filename, 'pdf') 
    clear fig_filename
jj = jj+1;
  end
clear jj
