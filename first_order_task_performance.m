%% first_order_task_performance
% This draws the plots in Figure 1B and 1D 
% (Figure 1C is not be plotted here as it requires the d-prime outputs from
% the analyses of local metacognitive efficiency, so it undertaken in that
% script)

% Get colour scheme
col = [[230, 10, 0]/255;... %Red
    [51, 140, 255]/255;... %Blue
    [60,179,113]/255;...%Green
    [255, 204, 204]/255;...%Light Red
    [100, 220, 255] / 255;...%Light Blue
    [152,251,152]/255;... %Light Green
    [85,107,47]/255];%Deep Green

% Set the parameters and choose the domain
  jj=1;
  while jj <5
x_var = age_single; % plots on age first
x_limits = [17.8 85];
x_ticks = age_groupmeans;
x_ticklabels = [{'18-27'},{'28-37'},{'38-47'},{'48-57'},{'58-67'},{'68+'}];

 if jj == 1
    fig_filename = 'Fig1Bi_diff_on_age';
    domain = 1; % 1 for mem, 2 for perc
    figure(11) 
    y_var = memory_variables.difflevel; 
    y_limits = [0,11];
 elseif jj == 2 
    fig_filename = 'Fig1Di_diffstd_on_age';
    domain = 1; 
    figure(15) 
    y_var = memory_variables.diffstd; 
    y_limits = [0,6];
 elseif jj == 3 
    fig_filename = 'Fig1Bii_diff_on_age';
    domain = 2; 
    figure(12) 
     y_var = perception_variables.difflevel; 
    y_limits = [0,15];
 elseif jj == 4 
    fig_filename = 'Fig1Dii_diffstd_on_age';
    domain = 2; 
    figure(16) 
    y_var = perception_variables.diffstd; 
    y_limits = [0,6];
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