%% Computes predicted probability of repetition during each trial using the DBM and correlates it with RT 

% Runs an exhaustive grid search of the param space for each subject and
% picks the parameter values that give best correlation with behavior. 

% Trials_Info{1}.TableDetails
%     Col_1: "Trial No. (** Not serial no.) A row with 0 means the trial was rejected based on eye tracking.)"
%     Col_2: "Cue Validity: 1=no_change, 2=valid, 3=invalid."
%     Col_3: "Cued Side: (-1)=Left cue, 1=Right Cue."
%     Col_4: "Angle of Change."
%     Col_5: "Change Flag for LEFT stimulus: 0=no_change, 1=change."
%     Col_6: "Change Flag for RIGHT stimulus: 0=no_change, 1=change."
%     Col_7: "Duration of Stimuli + Cue display (in seconds)."
%     Col_8: "Response Key: (-1)=Left Change, 0=No_change, 1=Right


%% Code
clear;
addpath('functions')
load ('data/tACS_40Hz_woETrej.mat')
side = ["Left","Right"];
sess = ["Sham","Stim","Post"];
%input type is only "choiceRA"
var_type=["choiceRA"];

rt_means=get_subject_means(leftPPC,rightPPC);
for iSide = 1:2
    fprintf("Side:%d\n",iSide);
    if(iSide == 1)
        Trials_Info = leftPPC.Trials_Info;
    else
        Trials_Info = rightPPC.Trials_Info;
    end
    
    for iSess = 1:3
        fprintf("Session:%d\n",iSess);
       
        for iSubject=1:length(Trials_Info)
            data_sub = Trials_Info{iSubject};
            fprintf("Subject:%d\n",iSubject);
            if(iSess==1)
                all_blocks = data_sub.Sham;
            elseif(iSess==2)
                all_blocks = data_sub.Stim;
            elseif(iSess==3)
                all_blocks = data_sub.Post;
            end
            
            [stats_choice, rt_all] = process_datablock(all_blocks,rt_means{iSubject});
            
            rvals_choice{iSide,iSess}(iSubject,:) = stats_choice.rlist;
            pvals_choice{iSide,iSess}(iSubject,:) = stats_choice.plist;
            alpha_choice{iSide,iSess}(iSubject) = stats_choice.best_params.alph;
            pmean_choice{iSide,iSess}(iSubject) = stats_choice.best_params.priormean;  
            pred_vec_choice{iSide,iSess}(iSubject,:) = stats_choice.pred_vec;
            pred_vec_choice_gridbest{iSide,iSess}(iSubject,:) = stats_choice.gridbest_pred_vec;
            rtvals{iSide,iSess}(iSubject,:) = rt_all;
            bin_vec_choice{iSide,iSess}(iSubject,:) = stats_choice.bin_vec;
           
            
        end
    end
end

save('results/compareGridbest_normDivideByMeanThenLog_'+var_type+'.mat')
disp("results saved")
%plot_DBM_gridbest


function [stats_choice, rt_all] = process_datablock(all_blocks, rt_mean)

% Unpack data from array

%% Look at various binary descriptors of current trial, and examine if forward predictions of these binary values influence RT.
%% In other words, is there evidence that subjects are in fact tracking and estimating likelihood that next trial is "valid", or "left_cued", or "change_trial".
for iBlock = 1:5
    fprintf("Block:%d\n",iBlock);
    ablock = all_blocks(:,:,iBlock);
    subj_choice=ablock(:,9);
    startTrial=2;
    for iTrial = startTrial:size(ablock,1)
        choice_rep(iTrial) = subj_choice(iTrial) == subj_choice(iTrial-1) & subj_choice(iTrial)~=5;
    end
    rt = ablock(startTrial:end,10)./rt_mean;
    for trl=1:length(rt)
        if rt(trl)~=0
            rt(trl)=log(rt(trl));
        end
    end
    choice_rep=choice_rep(startTrial:end);
     
    %% simple, un-fitted smoothing
    
    params=struct('alph', 0.7, 'prior_mean', 0.5, 'prior_scale', 20);

    block_idx = ((iBlock-1)*length(rt)+1:(iBlock)*length(rt));
    
    
    ptraj_choice(block_idx) = est_prob(choice_rep, params.alph, params.prior_mean, params.prior_scale);
    
    choice_rep_all(block_idx) = choice_rep;
 
    rt_all(block_idx) = rt;
    
    run_exhaustive = 1;
    if run_exhaustive == 1
        alph_list = [0.2:0.04:0.9 0.95 0.97 0.99];
        priormean_list = [0.05:0.04:0.95];
        pmscale_list = [20];
        
        for sc = 1:length(pmscale_list)
            for pp = 1:length(priormean_list)
                for aa=1:length(alph_list)
                    block_idx = ((iBlock-1)*length(rt)+1:(iBlock)*length(rt));
                    ptraj_choice_grid{aa,pp,sc}(block_idx) = est_prob(choice_rep, alph_list(aa), priormean_list(pp), pmscale_list(sc));
                end
            end
        end
    end
    
    
end

[stats_choice]  = est_and_corr(choice_rep_all, params, rt_all, ptraj_choice,ptraj_choice_grid);

end

function [stats] = est_and_corr(bin_vec, params, rt, ptraj, ptraj_grid)

stats.bin_vec=bin_vec;

[r, p] = corrcoef(bin_vec, rt);
stats.raw = struct('r', r(1,2), 'p', p(1,2));        %Correlation of binary input vector with rt

ptraj_pred = zeros(size(ptraj));
ptraj_pred(bin_vec==1) = ptraj(bin_vec==1);
ptraj_pred(bin_vec==0) = 1-ptraj(bin_vec==0);
% ptraj_pred = ptraj;
stats.pred_vec = ptraj;
[r, p] = corrcoef(ptraj_pred(rt~=0), rt(rt~=0));
stats.predictive = struct('r', r(1,2), 'p', p(1,2)); %Correlation of predictive input vector with rt <- value of forward prediction

%% SanityCheck: check that there is no target leakage; ptraj should ideally have no predictive power on actual label
% However, it is nonzero -- perhaps in high-frequency bit vectors with exact occupancy, there may still be some correlation.

[r, p] = corrcoef(bin_vec, ptraj');
stats.pred_vs_actual = struct('r', r(1,2), 'p', p(1,2));

%%% SanityCheck: Permute and smooth; this predictive estimate should have no
%%% correlation with RT.
permidx = randperm(length(bin_vec));
perm_bitvec = bin_vec(permidx); perm_rt = rt(permidx);
perm_ptraj = est_prob(perm_bitvec, params.alph, params.prior_mean, params.prior_scale);

[r, p] = corrcoef(perm_ptraj, perm_rt);
stats.permuted = struct('r', r(1,2), 'p', p(1,2));

%% combine both sanity checks: permuted and pred vs actual
[r, p] = corrcoef(perm_bitvec, perm_ptraj');
stats.permuted_predvsact = struct('r', r(1,2), 'p', p(1,2));

%% "best-effort" smoothing
%% Note: landscapes are weird, edgecases abound, params conflict with each other.
%% Quick check suggests there is no substantial improvement over randomly picked smoothing params; comment out for now.
run_exhaustive = 1;
if run_exhaustive == 1
    alph_list = [0.2:0.04:0.9 0.95 0.97 0.99];
    priormean_list = [0.05:0.04:0.95];
    pmscale_list = [20];
    
    for sc = 1:length(pmscale_list)
        for pp = 1:length(priormean_list)
            for aa=1:length(alph_list)
                ptraj_pred_grid{aa,pp,sc} = zeros(size(ptraj_grid{aa,pp,sc}));
                ptraj_pred_grid{aa,pp,sc}(bin_vec==1) = ptraj_grid{aa,pp,sc}(bin_vec==1);
                ptraj_pred_grid{aa,pp,sc}(bin_vec==0) = 1-ptraj_grid{aa,pp,sc}(bin_vec==0);
%                 ptraj_pred_grid{aa,pp,sc} = ptraj_grid{aa,pp,sc};
                [r, p] = corrcoef(ptraj_pred_grid{aa,pp,sc}(rt~=0), rt(rt~=0));
                rcache(aa,pp,sc) = r(1,2); pcache(aa,pp,sc) = p(1,2);
            end
        end
    end
    
    
    [~, idx] = min((rcache(:))); % max r (should correspond to most significant too).
    
    siz = size(rcache); [aidx, pidx, sidx] = ind2sub(siz, idx);
    
    alph = alph_list(aidx); % and p-value. (def < 0.05)
    pri = priormean_list(pidx);
    sca = pmscale_list(sidx);
    
    r = rcache(aidx,pidx,sidx); p = pcache(aidx,pidx,sidx);
    
    statsbest = struct('r', r, 'p', p);
    bestparams = struct('alph', alph, 'priormean', pri, 'scale', sca);
    stats.gridbest = statsbest;
    stats.best_params = bestparams;
    stats.gridbest_pred_vec = ptraj_grid{aidx,pidx,sidx};
else
    stats.names = {'raw' 'predictive' 'permuted' 'pred_vs_actual' 'permuted_predvsact'};
end

stats.names = {'raw' 'predictive', 'pred_vs_actual','gridbest'};

for i=1:length(stats.names)
    stats.rlist(i) = stats.(stats.names{i}).r;
    stats.plist(i) = stats.(stats.names{i}).p;
end

fprintf('.');

end


