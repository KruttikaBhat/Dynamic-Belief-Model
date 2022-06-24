%% Computes predicted probability of repetition during each trial using the DBM and correlates it with Accuracy


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
var_type = ["respEvent"];

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
            
            [stats_valid, stats_invalid, stats_nochange, stats_cue, acc_all] = process_datablock(all_blocks, var_type);
            
            rvals_valid{iSide,iSess}(iSubject,:) = stats_valid.rlist;
            rvals_invalid{iSide,iSess}(iSubject,:) = stats_invalid.rlist;
            rvals_nochange{iSide,iSess}(iSubject,:) = stats_nochange.rlist;
            rvals_cue{iSide,iSess}(iSubject,:) = stats_cue.rlist;
            
            pvals_valid{iSide,iSess}(iSubject,:) = stats_valid.plist;
            pvals_invalid{iSide,iSess}(iSubject,:) = stats_invalid.plist;
            pvals_nochange{iSide,iSess}(iSubject,:) = stats_nochange.plist;
            pvals_cue{iSide,iSess}(iSubject,:) = stats_cue.plist;
            
            alpha_valid{iSide,iSess}(iSubject) = stats_valid.best_params.alph;
            alpha_invalid{iSide,iSess}(iSubject) = stats_invalid.best_params.alph;
            alpha_nochange{iSide,iSess}(iSubject) = stats_nochange.best_params.alph;
            alpha_cue{iSide,iSess}(iSubject) = stats_cue.best_params.alph;     
            
            pmean_valid{iSide,iSess}(iSubject) = stats_valid.best_params.priormean;
            pmean_invalid{iSide,iSess}(iSubject) = stats_invalid.best_params.priormean;
            pmean_nochange{iSide,iSess}(iSubject) = stats_nochange.best_params.priormean;
            pmean_cue{iSide,iSess}(iSubject) = stats_cue.best_params.priormean;
            
            pred_vec_valid{iSide,iSess}(iSubject,:) = stats_valid.pred_vec;
            pred_vec_invalid{iSide,iSess}(iSubject,:) = stats_invalid.pred_vec;
            pred_vec_nochange{iSide,iSess}(iSubject,:) = stats_nochange.pred_vec;
            pred_vec_cue{iSide,iSess}(iSubject,:) = stats_cue.pred_vec;

            pred_vec_valid_gridbest{iSide,iSess}(iSubject,:) = stats_valid.gridbest_pred_vec;
            pred_vec_invalid_gridbest{iSide,iSess}(iSubject,:) = stats_invalid.gridbest_pred_vec;
            pred_vec_nochange_gridbest{iSide,iSess}(iSubject,:) = stats_nochange.gridbest_pred_vec;
            pred_vec_cue_gridbest{iSide,iSess}(iSubject,:) = stats_cue.gridbest_pred_vec;

            acc_vals{iSide,iSess}(iSubject,:) = acc_all;
            
            bin_vec_valid{iSide,iSess}(iSubject,:) = stats_valid.bin_vec;
            bin_vec_invalid{iSide,iSess}(iSubject,:) = stats_invalid.bin_vec;
            bin_vec_nochange{iSide,iSess}(iSubject,:) = stats_nochange.bin_vec;
            bin_vec_cue{iSide,iSess}(iSubject,:) = stats_cue.bin_vec;
        end
    end
end

save('results/Accuracy_compareGridbest_'+var_type+'.mat')
disp("results saved")
%plot_DBM_gridbest


function [stats_valid, stats_invalid, stats_nochange, stats_cue, acc_all] = process_datablock(all_blocks, var_type)

% Unpack data from array

%% Look at various binary descriptors of current trial, and examine if forward predictions of these binary values influence RT.
%% In other words, is there evidence that subjects are in fact tracking and estimating likelihood that next trial is "valid", or "left_cued", or "change_trial".
for iBlock = 1:5
    ablock = all_blocks(:,:,iBlock);
    fprintf("Block:%d\n",iBlock);
    if(contains(var_type,"actual"))
        valid_init = ablock(:, 2) == 2;
        invalid_init = ablock(:, 2) == 3;
        nochange_init = ablock(:, 2) ==1;
        cue_init = ablock(:,3) == -1;
       
    elseif(contains(var_type,"resp"))
        valid_init = ablock(:, 3) ==  ablock(:, 9);
        invalid_init = (ablock(:,3) == -1 & ablock(:,9)== 1) | (ablock(:,3) == 1 & ablock(:,9)== -1);
        nochange_init = ablock(:, 9) == 0;
        cue_init = ablock(:,3) == -1;        
    end
    if(var_type=="respRA1" || var_type=="respRA2")
        startTrial=2;
    else
        startTrial=1;
    end
    for iTrial = startTrial:size(ablock,1)
        
        %% For repetition tracking, uncomment this block
        %                 isvalid_cue_rep(iTrial) = isvalid_cue(iTrial) == isvalid_cue(iTrial-1);
        %                 isleft_cue_rep(iTrial) = isleft_cue(iTrial) == isleft_cue(iTrial-1);
        %                 ischange_trial_rep(iTrial) = ischange_trial(iTrial) == ischange_trial(iTrial-1);
        %% For Frequency tracking, uncomment this block
        %isvalid_cue_rep(iTrial) = isvalid_cue(iTrial);
        %isleft_cue_rep(iTrial) = isleft_cue(iTrial) ;
        %ischange_trial_rep(iTrial) = ischange_trial(iTrial) ;
        if(var_type=="respEvent")
            valid_rep(iTrial) = valid_init(iTrial);
            invalid_rep(iTrial) = invalid_init(iTrial);
            nochange_rep(iTrial) = nochange_init(iTrial);
            cue_rep(iTrial) = cue_init(iTrial);
        elseif(var_type=="respRA1")
            valid_rep(iTrial) = valid_init(iTrial) == valid_init(iTrial-1) & valid_init(iTrial)==1;
            invalid_rep(iTrial) = invalid_init(iTrial) == invalid_init(iTrial-1) & invalid_init(iTrial)==1;
            nochange_rep(iTrial) = nochange_init(iTrial) == nochange_init(iTrial-1) & nochange_init(iTrial)==1;
            cue_rep(iTrial) = ablock(iTrial,3)==ablock(iTrial-1,3);
        elseif(var_type=="respRA2")
            valid_rep(iTrial) = valid_init(iTrial) == valid_init(iTrial-1);
            invalid_rep(iTrial) = invalid_init(iTrial) == invalid_init(iTrial-1);
            nochange_rep(iTrial) = nochange_init(iTrial) == nochange_init(iTrial-1);
            cue_rep(iTrial) = ablock(iTrial,3)==ablock(iTrial-1,3);
        end
    end
    valid_rep=valid_rep(startTrial:end);
    invalid_rep=invalid_rep(startTrial:end);
    nochange_rep=nochange_rep(startTrial:end);
    cue_rep=cue_rep(startTrial:end);
    
    acc = (ablock(:,5)== 0 & ablock(:,6)== 0 & ablock(:,9) == 0) | (ablock(:,5)== 1 & ablock(:,9) == -1) | (ablock(:,6)== 1 & ablock(:,9) == 1);
    acc = double(acc);
    acc = acc(startTrial:end);
    %% simple, un-fitted smoothing
    
    params=struct('alph', 0.7, 'prior_mean', 0.5, 'prior_scale', 20);
    
    % estimate and return correlational statistics.
    % Check the correlation on RT for:
    %      raw binary vector,
    %      predictive prior expectation on trial, using the DBM,
    %      predictive prior expectation *after random permutation* (to confirm sequential effect)
    %      "best fit" smoothing by grid search over some param space for the DBM.
    
    block_idx = ((iBlock-1)*length(acc)+1:(iBlock)*length(acc));
    
    ptraj_valid(block_idx) = est_prob(valid_rep, params.alph, params.prior_mean, params.prior_scale);
    ptraj_invalid(block_idx) = est_prob(invalid_rep, params.alph, params.prior_mean, params.prior_scale);
    ptraj_nochange(block_idx) = est_prob(nochange_rep, params.alph, params.prior_mean, params.prior_scale);
    ptraj_cue(block_idx) = est_prob(cue_rep, params.alph, params.prior_mean, params.prior_scale);
    
    valid_rep_all(block_idx) = valid_rep;
    invalid_rep_all(block_idx) = invalid_rep;
    nochange_rep_all(block_idx) = nochange_rep;
    cue_rep_all(block_idx) = cue_rep;
    
    
    acc_all(block_idx) = acc;
    
    run_exhaustive = 1;
    if run_exhaustive == 1
        alph_list = [0.2:0.04:0.9 0.95 0.97 0.99];
        priormean_list = [0.05:0.04:0.95];
        pmscale_list = [20];
        
        for sc = 1:length(pmscale_list)
            for pp = 1:length(priormean_list)
                for aa=1:length(alph_list)
                    block_idx = ((iBlock-1)*length(acc)+1:(iBlock)*length(acc));
                    ptraj_valid_grid{aa,pp,sc}(block_idx) = est_prob(valid_rep, alph_list(aa), priormean_list(pp), pmscale_list(sc));
                    ptraj_invalid_grid{aa,pp,sc}(block_idx) = est_prob(invalid_rep, alph_list(aa), priormean_list(pp), pmscale_list(sc));
                    ptraj_nochange_grid{aa,pp,sc}(block_idx) = est_prob(nochange_rep, alph_list(aa), priormean_list(pp), pmscale_list(sc));
                    ptraj_cue_grid{aa,pp,sc}(block_idx) = est_prob(cue_rep, alph_list(aa), priormean_list(pp), pmscale_list(sc));

                end
            end
        end
    end
    
    
end

[stats_valid]  = est_and_corr(valid_rep_all, params, acc_all, ptraj_valid,ptraj_valid_grid);
[stats_invalid]  = est_and_corr(invalid_rep_all, params, acc_all, ptraj_invalid,ptraj_invalid_grid);
[stats_nochange]  = est_and_corr(nochange_rep_all, params, acc_all, ptraj_nochange,ptraj_nochange_grid);
[stats_cue]  = est_and_corr(cue_rep_all, params, acc_all, ptraj_cue,ptraj_cue_grid);
end

function [stats] = est_and_corr(bin_vec, params, acc, ptraj, ptraj_grid)

stats.bin_vec=bin_vec;

[r, p] = corrcoef(bin_vec, acc);
stats.raw = struct('r', r(1,2), 'p', p(1,2));        %Correlation of binary input vector with rt

ptraj_pred = zeros(size(ptraj));
ptraj_pred(bin_vec==1) = ptraj(bin_vec==1);
ptraj_pred(bin_vec==0) = 1-ptraj(bin_vec==0);
% ptraj_pred = ptraj;
stats.pred_vec = ptraj;
[r, p] = corrcoef(ptraj_pred, acc);
stats.predictive = struct('r', r(1,2), 'p', p(1,2)); %Correlation of predictive input vector with rt <- value of forward prediction

%% SanityCheck: check that there is no target leakage; ptraj should ideally have no predictive power on actual label
% However, it is nonzero -- perhaps in high-frequency bit vectors with exact occupancy, there may still be some correlation.

[r, p] = corrcoef(bin_vec, ptraj');
stats.pred_vs_actual = struct('r', r(1,2), 'p', p(1,2));

%%% SanityCheck: Permute and smooth; this predictive estimate should have no
%%% correlation with RT.
permidx = randperm(length(bin_vec));
perm_bitvec = bin_vec(permidx); perm_acc = acc(permidx);
perm_ptraj = est_prob(perm_bitvec, params.alph, params.prior_mean, params.prior_scale);

[r, p] = corrcoef(perm_ptraj, perm_acc);
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
                [r, p] = corrcoef(ptraj_pred_grid{aa,pp,sc}, acc);
                rcache(aa,pp,sc) = r(1,2); pcache(aa,pp,sc) = p(1,2);
            end
        end
    end
    
    
    [~, idx] = max((rcache(:))); % max r (should correspond to most significant too).
    
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


