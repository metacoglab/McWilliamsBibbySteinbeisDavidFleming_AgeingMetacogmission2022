%% preprocessing_and_exclusions
% This code processes the raw data from Metacogmission to produce the
% dataset for analysis
% Requires trials2counts

%% Get the raw data .mat files
% Load .mat files of the raw data set (n=304) for memory and perception
load ParticDemogs_and_globals % 
load Raw_data_memoryMetacogmission
load Raw_data_perceptionMetacogmission

for domain = 1:2; % The whole script runs twice to process the memory (1) and then perception (2) data 
if domain ==1 
    data_set = Raw_memory_data; 
else data_set = Raw_perception_data; 
end

numraw_trials =height(data_set); % Total raw number of TRIALS (totalling across all participants
numsubj_raw = max(data_set.cohort_ID);

% Add variable into table to indicate when first-order task was answered correctly 
data_set.answeredcorrectly = data_set.action==data_set.correctchoice;

% Delete some trials to account for burn in (trials 1 to burn_in_number) for each person 
burn_in_number =20; % enter number of trials for burn-in
clear rowsToDelete
rowsToDelete = data_set.rawtrialnumber< burn_in_number + 1;
%delete those trials from the raw data table
data_set(rowsToDelete,:) = [];
clear burn_in_number

%% Exclusions for long trials and long-trial subjects: 
% Trials with first order response time > 30s
longer30_first= data_set.actionTime>30000;
number_trials_longer30_first = sum(longer30_first); %Total number of trials to exclude
% Trials with confidence rating response times > 30s
longer30_conf= data_set.confidenceTime>30000;
number_trials_longer30_conf = sum(longer30_conf); %Total number of trials to exclude

% For info only, see if there are any trials where both conditions apply
clear number_trials_longer30_both
number_trials_longer30_both = sum(longer30_conf ==1 & longer30_first ==1); 

% Now find the row-numbers of the trials which need to be excluded
clear rowsToDelete
rowsToDelete = longer30_first==1 | longer30_conf==1;
clear longer30_first
clear longer30_conf

% delete those long trials from the raw data table
data_set(rowsToDelete,:) = [];
clear rowsToDelete

% Get rate of subject choosing correctly
subj_correct_rate=ones (numsubj_raw,1);
for j=1:numsubj_raw 
    subj_correct_rate(j)= mean(data_set.answeredcorrectly(data_set.cohort_ID==j)); 
 end
clear j

% Check that there are no under-performing subjects, getting less than 55%
% on the first-order task
 performance_min = 0.55; % this is the lowest acceptable first order accuracy rate
 first_performance_min=find(subj_correct_rate<performance_min);
number_subjs_to_exclude = numel(first_performance_min);
clear performance_min

clear numraw_trials

%% Re-number trials for each subject, now exclusions have been made
% Re-count number of trials, now exclusions have been made
numsubj= length(unique(data_set.cohort_ID)); % number of subjects after exclusion
% check that this is still 304
    if numsubj ~=numsubj_raw;
        error('some subjects need to be excluded')
    end
    clear numsubj_raw
num_trials=height(data_set); % total number of trials after exclusions

% Number of trials per subject
  for kk=1:numsubj
     trialcount_per_subj(kk)=sum(data_set.cohort_ID==kk);
 end
clear kk

% Re-number trials within each subject, now exlusions have been made
 % Do participant 1 first
 data_set.trialnumber(1:trialcount_per_subj(1))=1:trialcount_per_subj(1);
 clear next_subject_trials_begin_here
next_subject_trials_begin_here=trialcount_per_subj(1);
 % Now add the others in
  for kk=2:numsubj
        for j=1:trialcount_per_subj(kk) 
      data_set.trialnumber(next_subject_trials_begin_here + j)=j; %Number up the trials of single participant
        end
        clear j
      next_subject_trials_begin_here=next_subject_trials_begin_here + trialcount_per_subj(kk); 
            % increase the line number of where the last trial of the current subject was
  end
clear kk
 clear next_subject_trials_begin_here
 
%% Now add in performance variables 

%% Difficulty levels achieved by each subject on the staircase
% Mean difficulty
if domain ==1; % for memory
for kk=1:numsubj 
 difflevel(kk)= mean(data_set.difficulty(data_set.cohort_ID ==kk));
end
clear kk
else
% for perc, flip the number by subtracting from 15 so that higher values mean harder tasks:
for kk=1:numsubj 
 difflevel(kk)= 15.-mean(data_set.difficulty(data_set.cohort_ID ==kk));
end
clear kk
end
% ... and standard deviation of the staircased difficulty level for each
% participant
for kk=1:numsubj 
    diffstd(kk)=std(data_set.difficulty(data_set.cohort_ID ==kk));
end
clear kk

%% Calculate metacognitive bias (mean raw confidence per participant, irrespective of correctness
for kk=1:numsubj 
    bias(kk)=mean(data_set.confidence(data_set.cohort_ID ==kk));
end
clear kk

 output_variables =[]; % Put outputs into a structure
    output_variables.accuracy1storder = subj_correct_rate';
    output_variables.difflevel = difflevel;
    output_variables.diffstd = diffstd;
    output_variables.bias = bias;
    
    output_variables.numsubj = numsubj; 
    output_variables.num_trials = num_trials;
    output_variables.firstorder_underperformed = first_performance_min;
    output_variables.number_trials_longer30_first = number_trials_longer30_first;
    output_variables.number_trials_longer30_conf = number_trials_longer30_conf;
    output_variables.number_trials_longer30_both = number_trials_longer30_both;
    output_variables.number_subjs_to_exclude = number_subjs_to_exclude;
    output_variables.trialcount = trialcount_per_subj;
    
    if domain ==1;
        mem_output_variables = output_variables;
    else
        perc_output_variables = output_variables;
    end
%% Lastly, package confidence responses for use in analyses, by using trials2counts

% Variable for "response on left was the correct one" and "person chose
% left"
if domain == 1;
    correct_on_left = data_set.choiceleft == data_set.correctchoice;
    chose_left = data_set.choiceleft == data_set.action;
else
    correct_on_left = data_set.correctchoice =='A';
    chose_left = data_set.action =='A';
end

%% First make up the nR_S1 and nR_S2 vectors 
trials_nR_S1 = [];
trials_nR_S2 = [];

 %Bin the confidences
confidence_binned = [];
nbins=6; % number of bins for binned confidences
for kk=1:numsubj% 
confbins = quantile(data_set.confidence(data_set.cohort_ID==kk),nbins-1); % boundaries for a single subject for binning confidence into 6 groups
confbins = [-Inf confbins Inf];
for b = 1:length(confbins)-1
    % Change the raw confidences into the index of the bin they go in
    index = data_set.confidence(data_set.cohort_ID==kk) > confbins(b) & data_set.confidence(data_set.cohort_ID==kk) <= confbins(b+1);
    confidence_binned(index) = b;
end
clear b
conf = confidence_binned';% conf - vector of 1 x ntrials of confidence ratings taking values 1:nbins

% Now work out nR_S1 and nR_S2, putting subjects in on different lines, using trials2counts
[trials_nR_S1{kk} trials_nR_S2{kk}] = trials2counts(correct_on_left(data_set.cohort_ID==kk),chose_left(data_set.cohort_ID==kk),confidence_binned',nbins,0); 
clear confidence_binned
clear conf
clear index
end
clear kk 

    clear chose_left
    clear correct_on_left
    clear confbins
    clear nbins

% Re-name the variables
if domain == 1
    nR_S1_mem = trials_nR_S1;
    nR_S2_mem = trials_nR_S2;
    data_set_mem = data_set;  
    clear Raw_memory_data
else
    nR_S1_perc = trials_nR_S1;
    nR_S2_perc = trials_nR_S2;
    data_set_perc = data_set;
    clear Raw_perception_data
end

    clear trials_nR_S1
    clear trials_nR_S2

    clear subj_correct_rate;
    clear difflevel;
    clear diffstd;
    clear bias;
    
    clear numsubj; 
    clear num_trials;
    clear first_performance_min;
    clear number_trials_longer30_first;
    clear number_trials_longer30_conf;
    clear number_trials_longer30_both;
    clear number_subjs_to_exclude;
    clear trialcount_per_subj;
    clear output_variables
    
    clear data_set

end
clear domain
