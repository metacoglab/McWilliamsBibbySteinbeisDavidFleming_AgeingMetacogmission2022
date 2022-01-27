%% participant_characteristics
% Participant demographics analysis for McWilliams et al, 2022
% produces breakdown of participant characteristics shown in Table 1, and
% puts the means of each age group into a standalone vector for use in
% later plotting

load ParticDemogs_and_globals % get demographics data
partics=Partics_and_globals; % and rename it

age_single = partics.age_single';
age_group = partics.age_group';
Demogs_by_age_group = [];
for kk = 1:6
Demogs_by_age_group.groupsize(kk)= sum(age_group == kk);
Demogs_by_age_group.mean_age(kk)= mean(age_single(age_group == kk));
Demogs_by_age_group.std_age(kk)= std(age_single(age_group == kk));
Demogs_by_age_group.min_age(kk)= min(age_single(age_group == kk));
Demogs_by_age_group.max_age(kk)= max(age_single(age_group == kk));
Demogs_by_age_group.femaleness(kk) = sum(partics.gender(age_group==kk)=='female');
Demogs_by_age_group.maleness(kk) = sum(partics.gender(age_group==kk)=='male');
Demogs_by_age_group.prefernot(kk) = sum(partics.gender(age_group==kk)=='prefer_not_to_say');
end
clear kk

% Analyse sex balance for female v male within each age group

[tbl,chi2,p] = crosstab(Demogs_by_age_group.femaleness, Demogs_by_age_group.maleness);
Demogs_by_age_group.chi2stats.tbl = tbl;
Demogs_by_age_group.chi2stats.chi2 = chi2;
Demogs_by_age_group.chi2stats.p = p;

% Take mean ages for use in plots in other scripts
age_groupmeans = Demogs_by_age_group.mean_age;

clear tbl
clear chi2
clear p
clear partics
