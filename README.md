# McWilliamsBibbySteinbeisDavidFleming_AgeingMetacogmission2022

All analyses were done using Matlab and the scripts provided are for running in Matlab.

We provide the following 3 .mat files of anonymised behavioural data, to enable replication of analyses and rapid generation of the figures in the paper. 

1.	ParticDemogs_and_globals.mat contains participant age, sex and global confidence ratings for both the memory and perception tasks
2.	Raw_data_memoryMetacogmission contains the raw memory data from the metacognitive task (‘Metacogmission’)
3.	Raw_data_perceptionMetacogmission contains the raw perceptual data from the metacognitive task (‘Metacogmission’)

You may need to set the paths the scripts point to as the current directories. 
setwd("/xxx")
The hierarchical modelling requires JAGS, which seems to run best in version 3.4.0, as compatability issues can occur with later versions. More information about hierarchical estimation of metacognitive efficiency can be found at (https://github.com/metacoglab/HMeta-d) (Fleming, 2017).

The following scripts will undertake the analyses and need to be run in this sequence:

participant_characteristics
This takes the basic demographic data and analyses age and sex, as reported in the first section of the Results in the manuscript and at Table 1. 

preprocessing_and_exclusions 
This takes the 2 sets of metacognitive data (one for memory, one for perception), excludes bad trials (according to our pre-registered criteria) and returns and saves in your current directory the 2 processed data files. Additionally, it also derives simple metrics of performance and of metacognitive bias (mean confidence on rating trial-by-trial performance), before re-packaging the trial-by-trial confidence rating data for the hierchical Bayesian analyses undertaken in later scripts.

first_order_task_performance
This draws the plots in Figure 1B and 1D. (Figure 1C is not be plotted here as it requires the d-prime outputs from the analyses of local metacognitive efficiency, so it undertaken in that script)

local_metacognitive_efficiency_analysis
This runs the analyses of group local metacognitive efficiency, including hierarchical regression and hierarchical estimation of inter-domain correlation, presented in Table 2 and in Figure 2. 

bias_analysis
This produces the plots shown in Figure 3

globals_analysis
This produces the plots in Figure 4, as well as the regressions reported in the manuscript text

All the required helper function and text files are found at https://github.com/metacoglab/HMeta-d

License

This code is being released with a permissive open-source license. You should feel free to use or adapt the utility code as long as you follow the terms of the license, which are enumerated below. If you make use of or build on the behavioral analyses, we would appreciate that you cite the paper.

Copyright (c) 2021, Andrew McWilliams and Stephen Fleming

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
