%% globals_analysis
% This produces the plots in Figure 4, 
% as well as the regressions and t-tests reported in the manuscript text

% Get colour scheme
col = [[230, 10, 0]/255;... %Red
    [51, 140, 255]/255;... %Blue
    [60,179,113]/255;...%Green
    [255, 204, 204]/255;...%Light Red
    [100, 220, 255] / 255;...%Light Blue
    [152,251,152]/255;... %Light Green
    [85,107,47]/255];%Deep Green

% Set the parameters and choose the domain
    y_limits = [0,11]; % y axis to plot globals on
  jj=1;
  while jj <15
    if jj <5 % first the 4 plots with age groups indicated
    x_var = age_single; % plots on age first
    x_limits = [17.8 85];
    x_ticks = age_groupmeans;
    x_ticklabels = [{'18-27'},{'28-37'},{'38-47'},{'48-57'},{'58-67'},{'68+'}];
 if jj == 1
    fig_filename = 'Fig4Ai_mem_pre_on_age';
    glm_filename = 'glm_Fig4Ai';
    domain = 1; % 1 for mem, 2 for perc
    figure(41) 
    y_var = partics.PreMem; 
 elseif jj == 2 
    fig_filename = 'Fig4Aii_mem_post_on_age';
    glm_filename = 'glm_Fig4Aii';
    domain = 1; 
    figure(42) 
    y_var = partics.PostMem; 
 elseif jj == 3 
    fig_filename = 'Fig4Aiii_perc_pre_on_age';
    glm_filename = 'glm_Fig4Aiii';
    domain = 2; 
    figure(43) 
     y_var = partics.PrePerc; 
 elseif jj == 4 
    fig_filename = 'Fig4Aiv_perc_post_on_age';
    glm_filename = 'glm_Fig4Aiv';
    domain = 2; 
    figure(44) 
    y_var = partics.PostPerc; 
 end
end

if jj >4 % then the 4 plots of difficulty staircase level
 if jj == 5
    fig_filename = 'Fig4BAi_mem_pre_on_diff';
    glm_filename = 'glm_Fig4Bi';
    domain = 1; % 1 for mem, 2 for perc
    figure(45) 
    x_var = memory_variables.difflevel; % plots on staircase mean difficulty level
    x_limits = [0,12];
    y_var = partics.PreMem; 
 elseif jj == 6 
    fig_filename = 'Fig4Bii_mem_post_on_diff';
    glm_filename = 'glm_Fig4Bii';
    domain = 1; 
    figure(46) 
    x_var = memory_variables.difflevel; 
    x_limits = [0,12];
    y_var = partics.PostMem; 
 elseif jj == 7 
    fig_filename = 'Fig4Biii_perc_pre_on_diff';
    glm_filename = 'glm_Fig4Biii';
    domain = 2; 
    figure(47) 
    x_var = perception_variables.difflevel;
    x_limits = [0,15];
     y_var = partics.PrePerc; 
 elseif jj == 8 
    fig_filename = 'Fig4Biv_perc_post_on_diff';
    glm_filename = 'glm_Fig4Biv';
    domain = 2; 
    figure(48) 
    x_var = perception_variables.difflevel; 
    x_limits = [0,15];
    y_var = partics.PostPerc; 
 elseif jj == 9
    fig_filename = 'Globalnotplotted_mem_update_on_age';
    glm_filename = 'glm_Globalnotplotted_mem_update_on_age';
    domain = 1; % 1 for mem, 2 for perc
    figure(51) 
    x_var = age_single; 
    x_limits = [17.8 85];
    y_var = partics.PostMem-partics.PreMem; 
    y_limits = [-10,10];
 elseif jj == 10
    fig_filename = 'Globalnotplotted_perc_update_on_age';
    glm_filename = 'glm_Globalnotplotted_perc_update_on_age';
    domain = 2; 
    figure(52) 
    x_var = age_single; 
    x_limits = [17.8 85];
    y_var = partics.PostPerc-partics.PrePerc;    
    y_limits = [-10,10];
 elseif jj == 11
    fig_filename = 'Globalnotplotted_mem_update_on_diff';
    glm_filename = 'glm_Globalnotplotted_mem_update_on_diff';
    domain = 1; % 1 for mem, 2 for perc
    figure(53) 
    x_var = memory_variables.difflevel; % plots on staircase mean difficulty level
    x_limits = [0,12];
    y_var = partics.PostMem-partics.PreMem; 
    y_limits = [-10,10];
 elseif jj == 12
    fig_filename = 'Globalnotplotted_perc_update_on_diff';
    glm_filename = 'glm_Globalnotplotted_perc_update_on_diff';
    domain = 2; 
    figure(54) 
    x_var = perception_variables.difflevel; % plots on staircase mean difficulty level
    x_limits = [0,15];
    y_var = partics.PostPerc-partics.PrePerc; 
    y_limits = [-10,10];
 elseif jj == 13
    fig_filename = 'Globalnotplotted_percmean_on_memmean';
    glm_filename = 'glm_Globalnotplotted_percmean_on_memmean';
    domain = 3; % 1 for mem, 2 for perc
    figure(55) 
    x_var = (partics.PostMem+partics.PreMem)/2;
    x_limits = [0,12];
    y_var = (partics.PostPerc+partics.PrePerc)/2; 
    y_limits = [0,12];
 elseif jj == 14
    fig_filename = 'Globalnotplotted_update_perc_on_mem';
    glm_filename = 'glm_Globalnotplotted_update_perc_on_mem';
    domain = 3; 
    figure(56) 
    x_var = partics.PostMem-partics.PreMem; 
    x_limits = [-10,10];
    y_var = partics.PostPerc-partics.PrePerc; 
    y_limits = [-10,10];
 else
 end
end

% Draw figure
set(gcf, 'Position', [800 400 190 245],'Color',[1,1,1]);

box('off');
hold('all');

%% Scatter plot the data
scatter ( x_var, y_var, 12, 'MarkerEdgeColor', col(domain+3, :),'MarkerFaceColor', col(domain+3, :));
hold on
 
%% Add means for the 6 age groups
if jj<5
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
else
end
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
glm_outputs.b = b(:);
glm_outputs.sem = stats.se(:);
glm_outputs.p = stats.p(:);

save (glm_filename,'glm_outputs')
clear glm_filename
clear glm_outputs
clear b
clear dev
clear stats
clear b_unnorm
clear dev_unnorm
clear stats_unnorm

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

% Tigthen up margins
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
    
clear x_var
clear y_var

jj = jj+1;
  end
clear jj

% ttest comparing updates post-task to pre-task confidence within subjects
[memory_t.h,memory_t.p,memory_t.ci,memory_t.stats] = ttest (partics.PreMem, partics.PostMem);
save('memory_globalsupdate_ttest','memory_t')
[perception_t.h,perception_t.p,perception_t.ci,perception_t.stats] = ttest (partics.PrePerc, partics.PostPerc);
save('perception_globalsupdate_ttest','perception_t')