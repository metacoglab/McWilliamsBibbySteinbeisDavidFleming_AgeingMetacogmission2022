%% bias_analysis
% This produces the plots in Figure 3, 
% as well as the regression reported in the manuscript text

% Get colour scheme
col = [[230, 10, 0]/255;... %Red
    [51, 140, 255]/255;... %Blue
    [60,179,113]/255;...%Green
    [255, 204, 204]/255;...%Light Red
    [100, 220, 255] / 255;...%Light Blue
    [152,251,152]/255;... %Light Green
    [85,107,47]/255];%Deep Green

load ParticDemogs_and_globals
load outputs_memory
load outputs_perception
load age_means_by_group
age_single = Partics_and_globals.age_single;
age_group = Partics_and_globals.age_group;

% Set the parameters and choose the domain
jj =1;
while jj<3
    x_var = age_single; % plots on age first
    x_limits = [17.8 85];
    x_ticks = age_groupmeans;
    x_ticklabels = [{'18-27'},{'28-37'},{'38-47'},{'48-57'},{'58-67'},{'68+'}];
    y_limits = [-100,100];
    x_label ='age(years)';
    y_label ='metacognitive bias';
 if jj == 1
    fig_filename = 'Fig3Ai_mem_bias_on_age';
    glm_filename = 'glm_Fig3Ai';
    domain = 1; % 1 for mem, 2 for perc
    figure(31) 
    y_var = mem_output_variables.bias; 
    set(gcf, 'Position', [200 500 597 450],'Color',[1,1,1]);
 else 
    fig_filename = 'Fig3Aii_perc_bias_on_age';
    glm_filename = 'glm_Fig3Aii';
    domain = 2; 
    figure(32) 
    y_var = perc_output_variables.bias; 
    set(gcf, 'Position', [200 20 597 450],'Color',[1,1,1]);
end
% Draw figure

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

% GLM fitting and Add regression curve if wanted
cov = horzcat(normalize(x_var)); % normalise for the regression, though not for the plots
[b,dev,stats]= glmfit(cov,normalize(y_var), 'normal');

% Un-normalised version of glmfit for plots
[b_unnorm,dev_unnorm,stats_unnorm]= glmfit(x_var,y_var, 'normal');
hold on
polycoeffs = stats_unnorm.beta'; % coefficients of the fitted polynomial derived from glmfit
fittedcurve = polycoeffs(1) + polycoeffs(2).*x_limits; 
plot(x_limits, fittedcurve, 'k', 'LineWidth', 2);   
 
hold('all');

%% Saving
% glm_outputs.b = b(:);
% glm_outputs.sem = stats.se(:);
% glm_outputs.p = stats.p(:);

clear glm_filename
% clear glm_outputs

dim = [.8 .05 .9 .3]; % Now make text labels with betas and their p values
 if stats.p(2) >0.0001 
    annotation('textbox',dim,'String',sprintf('b = %.2f \np = %.2g',b(2),stats.p(2)) ,'FitBoxToText','on');
 else
    annotation('textbox',dim,'String',sprintf('b = %.2f \np < 0.0001',b(2)) ,'FitBoxToText','on');
 end

clear b
clear dev
clear stats
clear b_unnorm
clear dev_unnorm
clear stats_unnorm

xlabel(x_label)
ylabel(y_label)

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

% Tighten up margins
 tightInset = get(gca, 'TightInset');
position(1) = tightInset(1);
position(2) = tightInset(2);
position(3) = 1 - tightInset(1) - tightInset(3);
position(4) = 1 - tightInset(2) - tightInset(4);
set(axh, 'Position', position);

% savefig (gcf,fig_filename) %  now save the figures if wanted 
% saveas(gcf,fig_filename, 'pdf') % save pdf version if wanted
clear fig_filename

clear line_for_std
clear fittedcurve
clear tightInset  
clear polycoeffs
clear position
clear cov
clear x_var
clear y_var
% clear b_print
% clear p_print

jj = jj+1;
  end
clear jj
clear axh

%% Now plot perception bias on memory bias, for Figure 3B

 fig_filename = 'Fig3B_percbias_on_membias';
 glm_filename = 'glm_Fig3B';
 domain = 3; 
 figure(33) 
 x_limits = [-100,100];
 y_limits = [-100,100];
 x_var = mem_output_variables.bias;
 y_var = perc_output_variables.bias; 
    
% Draw figure
set(gcf, 'Position', [850 500 597 450],'Color',[1,1,1]);

box('off');
hold('all');

%% Scatter plot the data
scatter ( x_var, y_var, 12, 'MarkerEdgeColor', col(domain+3, :),'MarkerFaceColor', col(domain+3, :));
hold on
    
% GLM fitting to add regression lines
[b,dev,stats]= glmfit(normalize(x_var),normalize(y_var), 'normal');

% Un-normalised version of glmfit for plots
[b_unnorm,dev_unnorm,stats_unnorm]= glmfit(x_var,y_var, 'normal');
hold on
polycoeffs = stats_unnorm.beta'; % coefficients of the fitted polynomial derived from glmfit
fittedcurve = polycoeffs(1) + polycoeffs(2).*x_limits; 
plot(x_limits, fittedcurve, 'k', 'LineWidth', 2);   

dim = [.7 .05 .9 .3]; % Now make text labels with betas and their p values
 if stats.p(2) >0.0001 
    annotation('textbox',dim,'String',sprintf('b = %.2f \np = %.2g',b(2),stats.p(2)) ,'FitBoxToText','on');
 else
    annotation('textbox',dim,'String',sprintf('b = %.2f \np < 0.0001',b(2)) ,'FitBoxToText','on');
 end
 
hold('all');

%% Set axes
axh = gca;
axh.FontSize = 10;
axh.FontName = 'Arial';
% set(axh,'LooseInset',get(axh,'TightInset'));
axh.XAxisLocation = 'origin';
axh.YAxisLocation = 'origin';

x_label = 'memory';
y_label = 'perception';

xlim(x_limits)
ylim(y_limits)
xlabel(x_label)
ylabel(y_label)

%% Save glm outputs if you wish
% glm_outputs.b = b(:);
% glm_outputs.sem = stats.se(:);
% glm_outputs.p = stats.p(:);
% save (glm_filename,'glm_outputs')

clear glm_filename
clear glm_outputs
clear cov
clear b
clear dev
clear stats
clear b_unnorm
clear dev_unnorm
clear stats_unnorm
clear axh  
 % savefig (gcf,fig_filename) % now save the figure if required
 % saveas(gcf,fig_filename, 'pdf') 
clear fig_filename
    
clear x_var
clear y_var    
clear polycoeffs
clear fittedcurve
     
%% Finally, plot the regression coefficients for domain-generality of
% metacognitve bias within each of the 6 age groups

%%  6 age groups treated individually, with a regression of perception bias on memory bias within each
for kk = 1:6
    clear y_var
    clear cov
    clear b
    clear stats
    clear dev
    y_beta = normalize(perc_output_variables.bias(age_group==kk));
    cov = normalize(mem_output_variables.bias(age_group==kk));
    [b,dev,stats]= glmfit(cov, y_beta, 'normal');
    bias_DG_byAge_GLM_b(kk)= b(2);
    bias_DG_byAge_GLM_sem(kk)= stats.se(2);
    bias_DG_byAge_GLM_p(kk)= stats.p(2);
end

    fig_filename = 'Fig3C_biascorr_on_age';
    glm_filename = 'glm_Fig3C';
    domain = 3; 
    figure(34)
    set(gcf, 'Position', [850 20 597 450],'Color',[1,1,1]);
    x_limits = [17.8 85];
    x_ticks = age_groupmeans;
    x_ticklabels = [{'18-27'},{'28-37'},{'38-47'},{'48-57'},{'58-67'},{'68+'}];
    y_limits = [0.3,1];
%% Add s.d.s for the 6 age groups
for kk = 1:6
 clear line_for_std
line_for_std = line([age_groupmeans(kk),age_groupmeans(kk)],...
    [bias_DG_byAge_GLM_b(kk)- bias_DG_byAge_GLM_sem(kk),...
    bias_DG_byAge_GLM_b(kk)+ bias_DG_byAge_GLM_sem(kk)]);
    set(line_for_std,'Color',col(domain,:),'LineWidth',4);
hold on
end
clear kk
% Then betas themselves
plot (age_groupmeans, bias_DG_byAge_GLM_b, '-o',...
    'color', col(domain+4,:),...
    'LineWidth',4,...
    'MarkerSize',8,...
    'MarkerEdgeColor',col(domain+4,:),...
    'MarkerFaceColor',col(domain+4,:));
 hold on   

clear bias_DG_byAge_GLM_b
clear bias_DG_byAge_GLM_sem
clear bias_DG_byAge_GLM_p
 
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
x_label = 'age(years)';
y_label = 'regression betas';
xlabel(x_label)
ylabel(y_label)

%% Saving
% glm_outputs.b = bias_DG_byAge_GLM_b(:);
% glm_outputs.sem = bias_DG_byAge_GLM_sem(:);
% glm_outputs.p = bias_DG_byAge_GLM_p(:);
% save (glm_filename,'glm_outputs')
clear glm_filename
clear cov
clear dim
clear glm_outputs
clear b
clear dev
clear stats
clear b_unnorm
clear dev_unnorm
clear stats_unnorm
clear axh
clear line_for_std
 
%  savefig (gcf,fig_filename) % now save the figure 
%  saveas(gcf,fig_filename, 'pdf') 
clear fig_filename
    
clear x_var
clear y_var
clear x_label
clear x_limits
clear x_ticklabels
clear x_ticks
clear y_beta
clear y_group_means
clear y_group_std
clear y_label
clear y_limits
