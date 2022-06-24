ReadMe: Dynamic Belief Model (DBM) fit to Ankita's 2ADC (PPC-tACS) dataset, Sanjna's 2ADC (TMS) dataset and Varsha's dataset

(Last Updated: 24/06/2022)


Overview: 
The Dynamic Belief Model computes the predicted probability of an event based on the
sequence/history of events using a bayesian framework. This model can capture seemingly irrational
sequential effects in behavior, where subjects respond better/faster depending on the (uninformative)
recent sequence of stimuli. Here, we compute the model's predicted probability of a particular event 
on each trial and correlate that with their Reaction Times (RT) or Accuracy across trials. A significant
correlation might suggest that subjects are implicitly tracking that event, and that it influences 
their behavior. We can compare how tACS to PPC changes these history effects. 

Notes:
1)  Since the DBM is built for binary choice tasks, we have currently binarized the events/variables in the 2ADC. 
    The variables/events we use are: 	
	1. Event based logic:	
          a) valid: 1 if subject's response is valid, 0 otherwise
  	  b) invalid: 1 if subject's response is invalid, 0 otherwise
  	  c) nochange: 1 if subject's response is 'no change', 0 otherwise
          d) cue: 1 if left side was cued, 0 if right side was cued
	2. Repetition Alternation 1(AND) logic:
	  a) valid: 1 if subject's current response and previous response is valid, 0 otherwise
  	  b) invalid: 1 if subject's current response and previous response is invalid, 0 otherwise
  	  c) nochange: 1 if subject's current response and previous response is 'no change', 0 otherwise
          d) cue: 1 if left/right side was cued on current trial and previous trial, 0 if different sides were cued on current trial and previous trial
	3. Repetition Alternation 2(XOR) logic:
	  a) valid: 1 if subject's current response and previous response is valid or 'not valid' (i.e. can be invalid or nochange trial), 0 otherwise
  	  b) invalid: 1 if subject's current response and previous response is invalid or 'not invalid' (i.e valid or nochange trial)', 0 otherwise
  	  c) nochange: 1 if subject's current response and previous response is 'no change' or 'not no change' (i.e. valid or invalid trial), 0 otherwise
          d) cue: 1 if left/right side were cued on current trial and previous trial, 0 if different sides were cued on current trial and previous trial

2)  Model parameters:
        
        a) alpha: probability that p(t) = p(t-1) (where p is probability of event)
        b) prior distribution: mean & scale parameters of beta-distribution.
        

3)  The code takes in the sequence of these events, and computes the predicted probability of these events 
    occuring on each trial based on the model's assumptions and the parameter values. Then, it correlates 
    the pred. prob. with RT/accuracy across trials. 
    
4)  Note that we currently do not estimate the model parameters. We do a gridsearch of the parameter values 
    that give the best correlations of predicted probabilities with RT/accuracy for each subject. 

A-Code file present in Ankita's subfolder
S-Code file present in Sanjna's subfolder
V-Code file present in Varsha's subfolder
C-Code file present in Combined datasets subfolder - this has results that are pooled across the 3 datasets.

Code:

1. corr_pred_delta_sens_bias_delta (ASVC): Gets the DBM 'r' (correlation between DBM predicted probability and RT) for the valid and invalid trials separately and computes delta DBM 'r'. This is plotted against delta d'/cc.
2. corr_pred_sens_bias (ASVC): Plots the DBM 'r' against the d'/cc.
3. corr_pred_sens_bias_delta (ASVC): Plots the DBM 'r' against delta d'/cc.
4. corr_pred_sens_bias_valInvAvg (C): Plots the DBM 'r' against d' in which d' valid and d' invalid are averaged.
5. corr_val_inv_DBMr (ASV): Plots the DBM 'r' of the valid variable against the DBM 'r' of the invalid variable.
6. DBM_gridbest_accuracy (A): Applies DBM model subject-session-blockwise to the 4 variables based on the input encoding that is specified. Correlation is taken between the predicted probabilities and accuracy. Gridbest method is applied.
7. DBM_gridbest_RT (ASV): Applies DBM model subject-session-blockwise to the 4 variables based on the input encoding that is specified. Correlation is taken between the predicted probabilities and RT. Gridbest method is applied.
8. dbm_pred_seq_effects (A): Bins trials based on recent sequence of events and plots the mean RT/1-DBM predicted probability(for RT comparison)/DBM predicted probability(for accuracy comparison) for each bin. (Based on figure 1A (Yu,Cohen 2008)) 
9. get_data (SV): Gets the required data from the subject's raw behavioral data folders. (all_blocks_data.mat)
10. get_subject_means (ASV): Function for getting the subject means, in order to do a normalisation(divide by subject mean-then-log) in the DBM code.
11. plot_DBM_gridbest (A): Plot the DBM 'r' in ascending order (subject number in x-axis does not reflect the actual subject number).
12. RT_pred_effects (A): Bins trials based on recent sequence of events and plots the mean RT for each bin. (Based on figure 1A (Yu,Cohen 2008)). Same result as dbm_pred_seq_effects if 'RT' is set as the beh_metric.

Combined choices DBM model codes: DBM is fit for subject's responses. 1 is if the subject's current and previous response is the same, 0 otherwise.
1. corr_pred_sens_bias_combinedChoices (A): Plots the DBM 'r' against the d'/cc for the combined choice DBM model.
2. DBM_gridbest_RT_combinedChoices (A):  Applies DBM model subject-session-blockwise to the choice RA variable. Correlation is taken between the predicted probabilities and RT. Gridbest method is applied.
3. dbm_pred_seq_effects_combinedChoices (A): Bins trials based on recent sequence of events and plots the mean RT/1-DBM predicted probability(for RT comparison) for each bin. (Based on figure 1A (Yu,Cohen 2008)) 
4. plot_DBM_gridbest_combinedChoices (A): Plot the DBM 'r' in ascending order (subject number in x-axis does not reflect the actual subject number).


Additional Notes:
1) Initial code sent by Dr. Pradeep Shenoy (see folder "PradeepShenoy_DBM_code" for the code he sent)
2) /data
    i) Ankita : raw data: tACS_40Hz_woETrej.mat, d'/bias: Indiv_Extracted_Data_(Left/Right)_ettag_noET.mat. I have also put the d' and cc in the .mat file called Ankita_perf_metrics.mat for consistency.
    ii) Sanjna : raw_data: all_blocks_data.mat, d'/cc: Sanjna_perf_metrics.mat
    iii) Varsha: 
	a) Sham-1(odd blocks),Sham-2(even blocks) : raw_data: all_blocks_data.mat, d'/cc: Varsha_oddeven_perf_metrics.mat
	b) Sham(all 10 blocks are called as Sham session) : raw_data: all_blocks_data_10blocks.mat, d'/cc: Varsha_10Blocks_perf_metrics.mat
3) perf_metrics file was obtained with the code 'get_perf_metrics'. This reads the json file specified. The json file is obtained with the python code 'get_Ankita_data.py','get_Sanjna_data.py', 'get_Varsha_data.py'.

References:
Sequential effects: Superstition or rational behavior? - Yu & Cohen, NeurIPS 2008
