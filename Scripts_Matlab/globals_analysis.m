%% globals_analysis
% This produces the plots in Figure 4, 
% as well as the regressions and t-tests reported in the manuscript text
% The final plots (Figures 51-6) are not drawn in the manusctipt, but are reported only
% in the text, but are left here for interest

load ParticDemogs_and_globals % get demographics data
load age_means_by_group
load outputs_memory
load outputs_perception

partics=Partics_and_globals; % and rename it
age_single = partics.age_single;
age_group = partics.age_group;

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
    y_label = 'confidence';
  jj=1;
  while jj <15
    if jj <5 % first the 4 plots with age groups indicated
    x_var = age_single; % plots on age first
    x_limits = [17.8 85];
    x_ticks = age_groupmeans;
    x_ticklabels = [{'18-27'},{'28-37'},{'38-47'},{'48-57'},{'58-67'},{'68+'}];
    x_label = 'age group (years)'; 
 if jj == 1
    fig_filename = 'Fig4Ai_mem_pre_on_age'; % File names if you with to save figures
    glm_filename = 'glm_Fig4Ai'; % File names if you with to save glm outputs
    domain = 1; % 1 for mem, 2 for perc
    figure(41) 
    y_var = partics.PreMem; 
    % Draw figure
set(gcf, 'Position', [100 400 240 245],'Color',[1,1,1]);
 elseif jj == 2 
    fig_filename = 'Fig4Aii_mem_post_on_age';
    glm_filename = 'glm_Fig4Aii';
    domain = 1; 
    figure(42) 
    y_var = partics.PostMem; 
    % Draw figure
  set(gcf, 'Position', [350 400 240 245],'Color',[1,1,1]);
 elseif jj == 3 
    fig_filename = 'Fig4Aiii_perc_pre_on_age';
    glm_filename = 'glm_Fig4Aiii';
    domain = 2; 
    figure(43) 
     y_var = partics.PrePerc; 
  set(gcf, 'Position', [100 100 240 245],'Color',[1,1,1]);
 elseif jj == 4 
    fig_filename = 'Fig4Aiv_perc_post_on_age';
    glm_filename = 'glm_Fig4Aiv';
    domain = 2; 
    figure(44) 
    y_var = partics.PostPerc; 
  set(gcf, 'Position', [350 100 240 245],'Color',[1,1,1]);
 end

    elseif jj >4 % then the 4 plots of difficulty staircase level
    x_label = 'difficulty level (a.u.)';
 if jj == 5
    fig_filename = 'Fig4BAi_mem_pre_on_diff';
    glm_filename = 'glm_Fig4Bi';
    domain = 1; % 1 for mem, 2 for perc
    figure(45) 
    x_var = mem_output_variables.difflevel; % plots on staircase mean difficulty level
    x_limits = [0,12];
    x_ticks = 0:2:12;
    x_ticklabels = x_ticks;
    y_var = partics.PreMem; 
    set(gcf, 'Position', [600 400 240 245],'Color',[1,1,1]);
 elseif jj == 6 
    fig_filename = 'Fig4Bii_mem_post_on_diff';
    glm_filename = 'glm_Fig4Bii';
    domain = 1; 
    figure(46) 
    x_var = mem_output_variables.difflevel; 
    x_limits = [0,12];
    x_ticks = 0:2:12;
    x_ticklabels = x_ticks;
    y_var = partics.PostMem; 
    set(gcf, 'Position', [600 100 240 245],'Color',[1,1,1]);
 elseif jj == 7 
    fig_filename = 'Fig4Biii_perc_pre_on_diff';
    glm_filename = 'glm_Fig4Biii';
    domain = 2; 
    figure(47) 
    x_var = perc_output_variables.difflevel;
    x_limits = [0,15];
    x_ticks= [0;2;4;6;8;10;12;14];
    x_ticklabels = x_ticks;
    y_var = partics.PrePerc; 
    set(gcf, 'Position', [850 400 240 245],'Color',[1,1,1]);
 elseif jj == 8 
    fig_filename = 'Fig4Biv_perc_post_on_diff';
    glm_filename = 'glm_Fig4Biv';
    domain = 2; 
    figure(48) 
    x_var = perc_output_variables.difflevel; 
    x_limits = [0,15];
    x_ticks= [0;2;4;6;8;10;12;14];
    x_ticklabels = x_ticks;
    y_var = partics.PostPerc; 
    set(gcf, 'Position', [850 100 240 245],'Color',[1,1,1]);
 elseif jj == 9
    fig_filename = 'Globalnotplotted_mem_update_on_age';
    glm_filename = 'glm_Globalnotplotted_mem_update_on_age';
    domain = 1; % 1 for mem, 2 for perc
    figure(51) 
    x_var = age_single;
    x_label = 'age (years)'; 
    x_limits = [17.8 85];
    x_ticks = 20:10:80;
    x_ticklabels = x_ticks;
    y_var = partics.PostMem-partics.PreMem; 
    y_limits = [-10,10];
    set(gcf, 'Position', [1400 200 240 245],'Color',[1,1,1]);
 elseif jj == 10
    fig_filename = 'Globalnotplotted_perc_update_on_age';
    glm_filename = 'glm_Globalnotplotted_perc_update_on_age';
    domain = 2; 
    figure(52) 
    x_var = age_single; 
    x_label = 'age (years)'; 
    x_limits = [17.8 85];
    x_ticks = 20:10:80;
    x_ticklabels = x_ticks;
    y_var = partics.PostPerc-partics.PrePerc;    
    y_limits = [-10,10];
    set(gcf, 'Position', [1420 200 240 245],'Color',[1,1,1]);
 elseif jj == 11
    fig_filename = 'Globalnotplotted_mem_update_on_diff';
    glm_filename = 'glm_Globalnotplotted_mem_update_on_diff';
    domain = 1; % 1 for mem, 2 for perc
    figure(53) 
    x_var = mem_output_variables.difflevel; % plots on staircase mean difficulty level
    x_label = 'difficulty level (a.u.)';
    x_limits = [0,12];
    x_ticks = 0:2:12;
    x_ticklabels = x_ticks;
    y_var = partics.PostMem-partics.PreMem; 
    y_limits = [-10,10];
    set(gcf, 'Position', [1440 200 240 245],'Color',[1,1,1]);
 elseif jj == 12
    fig_filename = 'Globalnotplotted_perc_update_on_diff';
    glm_filename = 'glm_Globalnotplotted_perc_update_on_diff';
    domain = 2; 
    figure(54) 
    x_var = perc_output_variables.difflevel; % plots on staircase mean difficulty level
    x_label = 'difficulty level (a.u.)';
    x_limits = [0,15];
    x_ticks= [0;2;4;6;8;10;12;14];
    x_ticklabels = x_ticks;
    y_var = partics.PostPerc-partics.PrePerc; 
    y_limits = [-10,10];
    set(gcf, 'Position', [1460 200 240 245],'Color',[1,1,1]);
 elseif jj == 13
    fig_filename = 'Globalnotplotted_percmean_on_memmean';
    glm_filename = 'glm_Globalnotplotted_percmean_on_memmean';
    domain = 3; % 1 for mem, 2 for perc
    figure(55) 
    x_var = (partics.PostMem+partics.PreMem)/2;
    x_label = 'memory confidence';
    x_limits = [0,12];
    x_ticks = 0:2:12;
    x_ticklabels = x_ticks;
    y_var = (partics.PostPerc+partics.PrePerc)/2; 
    y_label = 'perception confidence';
    y_limits = [0,12];
    set(gcf, 'Position', [1480 200 240 245],'Color',[1,1,1]);
 elseif jj == 14
    fig_filename = 'Globalnotplotted_update_perc_on_mem';
    glm_filename = 'glm_Globalnotplotted_update_perc_on_mem';
    domain = 3; 
    figure(56) 
    x_var = partics.PostMem-partics.PreMem; 
    x_label = 'memory confidence';
    x_limits = [-10,10];
    x_ticks = -10:5:10;
    x_ticklabels = x_ticks;
    y_var = partics.PostPerc-partics.PrePerc; 
    y_label = 'perception confidence';
    y_limits = [-10,10];
    set(gcf, 'Position', [1500 200 240 245],'Color',[1,1,1]);
    % Draw figure
 else
 end
end

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

% glm_outputs.b = b(:); % Option for use if saving
% glm_outputs.sem = stats.se(:);
% glm_outputs.p = stats.p(:);

% Now make text labels with betas and their p values
if ismember(jj,[1:4])==1 % adjust loction on y-axis slightly according to the plot
    dim = [.2 .13 .9 .3]; 
elseif ismember(jj,[9,10])==1
    dim = [.2 .08 .9 .3]; 
else
    dim = [.2 .05 .9 .3]; 
end

 if stats.p(2) >0.0001 
    annotation('textbox',dim,'String',sprintf('b = %.2f \np = %.2g',b(2),stats.p(2)) ,'FitBoxToText','on');
 else
    annotation('textbox',dim,'String',sprintf('b = %.2f \np < 0.0001',b(2)) ,'FitBoxToText','on');
 end

% save (glm_filename,'glm_outputs')
% clear glm_filename
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
%set(axh,'LooseInset',get(axh,'TightInset'));

% Angle the tick labels when age groups are on the tick labels
if ismember(jj,[1:4,9,10])==1
xtickangle(45)
else
end

xlim(x_limits)
ylim(y_limits)
xlabel(x_label)
ylabel(y_label)

% Tigthen up margins
 tightInset = get(gca, 'TightInset');
position(1) = tightInset(1);
position(2) = tightInset(2);
position(3) = 1 - tightInset(1) - tightInset(3);
position(4) = 1 - tightInset(2) - tightInset(4);
set(axh, 'Position', position);

%    savefig (gcf,fig_filename) % now save the figure if desired
%    saveas(gcf,fig_filename, 'pdf') 
    clear fig_filename
    
clear x_var
clear y_var

jj = jj+1;
  end
clear jj

clear axh
clear cov
clear dim
clear domain
clear fittedcurve
clear glm_filename
clear line_for_std
clear polycoeffs
clear tightInset
clear x_label
clear x_limits
clear x_ticklabels
clear x_ticks
clear y_group_std
clear y_label
clear y_limits
clear y_group_means

%% Finally, the t-tests reported in the paper 
% ttest comparing updates post-task to pre-task confidence within subjects
[memory_t.h,memory_t.p,memory_t.ci,memory_t.stats] = ttest (partics.PreMem, partics.PostMem);
% save('memory_globalsupdate_ttest','memory_t')
[perception_t.h,perception_t.p,perception_t.ci,perception_t.stats] = ttest (partics.PrePerc, partics.PostPerc);
% save('perception_globalsupdate_ttest','perception_t')
clear memory_t
clear perception_t
