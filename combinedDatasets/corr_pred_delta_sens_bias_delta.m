
%% code
clear;
type="respEvent";
sess=["Sham"];
plot_savepath="results/corr_pred-delta_";

%get Ankita's data
load("../DBM_fit_Ankita_data/results/compareGridbest_"+type+".mat") % Load the DBM results file to be plotted
load("../DBM_fit_Ankita_data/data/Ankita_perf_metrics.mat")
load("../DBM_fit_Ankita_data/data/tACS_40Hz_woETrej.mat")

dataset=1;
for iSide = 1:2
    if iSide==1
        trial_data=leftPPC.Trials_Info;
    else
        trial_data=rightPPC.Trials_Info;
    end
    iSess=1;
        
    for iSubject = 1:26

        sensDelta{dataset}((iSide-1)*26+iSubject) = mean(d_val(iSubject,(iSide-1)*3+iSess,:))-mean(d_inv(iSubject,(iSide-1)*3+iSess,:));
        cc_biasDelta{dataset}((iSide-1)*26+iSubject) = cc_val(iSubject,(iSide-1)*3+iSess)-cc_inv(iSubject,(iSide-1)*3+iSess);

        rt=rtvals{iSide,iSess}(iSubject,:);

        if iSess==1
            sub_data=trial_data{iSubject}.Sham;
        elseif iSess==2
            sub_data=trial_data{iSubject}.Stim;
        elseif iSess==3
            sub_data=trial_data{iSubject}.Post;
        end
        valid_init=[];
        invalid_init=[];
        if(type=="respRA1" || var_type=="respRA2")
            startTrial=2;
        else
            startTrial=1;
        end
        for block=1:size(sub_data,3)        
            ablock=sub_data(:,:,block);
            valid_init = [valid_init (ablock(startTrial:end, 2) == 2)'];
            invalid_init = [invalid_init (ablock(startTrial:end, 2) == 3)'];     
        end
        pred_vec_valid=pred_vec_valid_gridbest{iSide,iSess}(iSubject,:);
        pred_vec_invalid=pred_vec_invalid_gridbest{iSide,iSess}(iSubject,:);
        pred_vec_nochange=pred_vec_nochange_gridbest{iSide,iSess}(iSubject,:);
        pred_vec_cue=pred_vec_cue_gridbest{iSide,iSess}(iSubject,:);

        vec_valid=bin_vec_valid{iSide,iSess}(iSubject,:);
        vec_invalid=bin_vec_invalid{iSide,iSess}(iSubject,:);
        vec_nochange=bin_vec_nochange{iSide,iSess}(iSubject,:);
        vec_cue=bin_vec_cue{iSide,iSess}(iSubject,:);

        ptraj_pred_valid = zeros(size(pred_vec_valid));
        ptraj_pred_valid(vec_valid==1) = pred_vec_valid(vec_valid==1);
        ptraj_pred_valid(vec_valid==0) = 1-pred_vec_valid(vec_valid==0); 
        [r_val, p_val] = corrcoef(ptraj_pred_valid(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
        [r_inv, p_inv] = corrcoef(ptraj_pred_valid(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
        predCorr_val{dataset}((iSide-1)*26+iSubject) = r_val(1,2)-r_inv(1,2);

        ptraj_pred_invalid = zeros(size(pred_vec_invalid));
        ptraj_pred_invalid(vec_invalid==1) = pred_vec_invalid(vec_invalid==1);
        ptraj_pred_invalid(vec_invalid==0) = 1-pred_vec_invalid(vec_invalid==0);           
        [r_val, p_val] = corrcoef(ptraj_pred_invalid(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
        [r_inv, p_inv] = corrcoef(ptraj_pred_invalid(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
        predCorr_inval{dataset}((iSide-1)*26+iSubject) = r_val(1,2)-r_inv(1,2);

        ptraj_pred_nochange = zeros(size(pred_vec_nochange));
        ptraj_pred_nochange(vec_nochange==1) = pred_vec_nochange(vec_nochange==1);
        ptraj_pred_nochange(vec_nochange==0) = 1-pred_vec_nochange(vec_nochange==0);           
        [r_val, p_val] = corrcoef(ptraj_pred_nochange(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
        [r_inv, p_inv] = corrcoef(ptraj_pred_nochange(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
        predCorr_nochange{dataset}((iSide-1)*26+iSubject) = r_val(1,2)-r_inv(1,2);

        ptraj_pred_cue = zeros(size(pred_vec_cue));
        ptraj_pred_cue(vec_cue==1) = pred_vec_cue(vec_cue==1);
        ptraj_pred_cue(vec_cue==0) = 1-pred_vec_cue(vec_cue==0);           
        [r_val, p_val] = corrcoef(ptraj_pred_cue(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
        [r_inv, p_inv] = corrcoef(ptraj_pred_cue(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
        predCorr_cue{dataset}((iSide-1)*26+iSubject) = r_val(1,2)-r_inv(1,2);


    end

end   

%get Sanjna's data
load("../DBM_fit_Sanjna_data/results/compareGridbest_"+type+".mat") % Load the DBM results file to be plotted
load("../DBM_fit_Sanjna_data/data/Sanjna_perf_metrics.mat")
load("../DBM_fit_Sanjna_data/data/all_blocks_data.mat")

iSess=1;
dataset=2;
for iSubject = 1:28

    sensDelta{dataset}(iSubject) = mean(d_val(iSubject,iSess,:))-mean(d_inv(iSubject,iSess,:));
    cc_biasDelta{dataset}(iSubject) = cc_val(iSubject,iSess)-cc_inv(iSubject,iSess);

    rt=rtvals{iSess}(iSubject,:);
    sub_data = Trials_info{iSubject,iSess};
    valid_init=[];
    invalid_init=[];
    if(type=="respRA1" || var_type=="respRA2")
        startTrial=2;
    else
        startTrial=1;
    end
    for block=1:size(sub_data,3)    
        ablock=sub_data(:,:,block);
        valid_init = [valid_init (ablock(startTrial:end, 3) == 2)'];
        invalid_init = [invalid_init (ablock(startTrial:end, 3) == 3)'];     
    end
    pred_vec_valid=pred_vec_valid_gridbest{iSess}(iSubject,:);
    pred_vec_invalid=pred_vec_invalid_gridbest{iSess}(iSubject,:);
    pred_vec_nochange=pred_vec_nochange_gridbest{iSess}(iSubject,:);
    pred_vec_cue=pred_vec_cue_gridbest{iSess}(iSubject,:);

    vec_valid=bin_vec_valid{iSess}(iSubject,:);
    vec_invalid=bin_vec_invalid{iSess}(iSubject,:);
    vec_nochange=bin_vec_nochange{iSess}(iSubject,:);
    vec_cue=bin_vec_cue{iSess}(iSubject,:);

    ptraj_pred_valid = zeros(size(pred_vec_valid));
    ptraj_pred_valid(vec_valid==1) = pred_vec_valid(vec_valid==1);
    ptraj_pred_valid(vec_valid==0) = 1-pred_vec_valid(vec_valid==0); 
    [r_val, p_val] = corrcoef(ptraj_pred_valid(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
    [r_inv, p_inv] = corrcoef(ptraj_pred_valid(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
    predCorr_val{dataset}(iSubject) = r_val(1,2)-r_inv(1,2);

    ptraj_pred_invalid = zeros(size(pred_vec_invalid));
    ptraj_pred_invalid(vec_invalid==1) = pred_vec_invalid(vec_invalid==1);
    ptraj_pred_invalid(vec_invalid==0) = 1-pred_vec_invalid(vec_invalid==0);           
    [r_val, p_val] = corrcoef(ptraj_pred_invalid(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
    [r_inv, p_inv] = corrcoef(ptraj_pred_invalid(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
    predCorr_inval{dataset}(iSubject) = r_val(1,2)-r_inv(1,2);

    ptraj_pred_nochange = zeros(size(pred_vec_nochange));
    ptraj_pred_nochange(vec_nochange==1) = pred_vec_nochange(vec_nochange==1);
    ptraj_pred_nochange(vec_nochange==0) = 1-pred_vec_nochange(vec_nochange==0);           
    [r_val, p_val] = corrcoef(ptraj_pred_nochange(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
    [r_inv, p_inv] = corrcoef(ptraj_pred_nochange(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
    predCorr_nochange{dataset}(iSubject) = r_val(1,2)-r_inv(1,2);

    ptraj_pred_cue = zeros(size(pred_vec_cue));
    ptraj_pred_cue(vec_cue==1) = pred_vec_cue(vec_cue==1);
    ptraj_pred_cue(vec_cue==0) = 1-pred_vec_cue(vec_cue==0);           
    [r_val, p_val] = corrcoef(ptraj_pred_cue(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
    [r_inv, p_inv] = corrcoef(ptraj_pred_cue(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
    predCorr_cue{dataset}(iSubject) = r_val(1,2)-r_inv(1,2);

end


%get Varsha's data
load("../DBM_fit_Varsha_data/results/10Blocks/compareGridbest_"+type+".mat") % Load the DBM results file to be plotted
load("../DBM_fit_Varsha_data/data/Varsha_10Blocks_perf_metrics.mat")
load ('../DBM_fit_Varsha_data/data/all_blocks_data_10blocks.mat')

iSess=1;
dataset=3;
for iSubject = 1:22

    sensDelta{dataset}(iSubject) = mean(d_val(iSubject,iSess,:))-mean(d_inv(iSubject,iSess,:));
    cc_biasDelta{dataset}(iSubject) = cc_val(iSubject,iSess)-cc_inv(iSubject,iSess);

    rt=rtvals{iSess}(iSubject,:);
    sub_data = Trials_info{iSubject,iSess};
    valid_init=[];
    invalid_init=[];
    if(type=="respRA1" || var_type=="respRA2")
        startTrial=2;
    else
        startTrial=1;
    end
    for block=1:size(sub_data,3)   
        ablock=sub_data(:,:,block);
        valid_init = [valid_init (ablock(startTrial:end, 3) == 2)'];
        invalid_init = [invalid_init (ablock(startTrial:end, 3) == 3)'];     
    end
    pred_vec_valid=pred_vec_valid_gridbest{iSess}(iSubject,:);
    pred_vec_invalid=pred_vec_invalid_gridbest{iSess}(iSubject,:);
    pred_vec_nochange=pred_vec_nochange_gridbest{iSess}(iSubject,:);
    pred_vec_cue=pred_vec_cue_gridbest{iSess}(iSubject,:);

    vec_valid=bin_vec_valid{iSess}(iSubject,:);
    vec_invalid=bin_vec_invalid{iSess}(iSubject,:);
    vec_nochange=bin_vec_nochange{iSess}(iSubject,:);
    vec_cue=bin_vec_cue{iSess}(iSubject,:);

    ptraj_pred_valid = zeros(size(pred_vec_valid));
    ptraj_pred_valid(vec_valid==1) = pred_vec_valid(vec_valid==1);
    ptraj_pred_valid(vec_valid==0) = 1-pred_vec_valid(vec_valid==0); 
    [r_val, p_val] = corrcoef(ptraj_pred_valid(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
    [r_inv, p_inv] = corrcoef(ptraj_pred_valid(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
    predCorr_val{dataset}(iSubject) = r_val(1,2)-r_inv(1,2);

    ptraj_pred_invalid = zeros(size(pred_vec_invalid));
    ptraj_pred_invalid(vec_invalid==1) = pred_vec_invalid(vec_invalid==1);
    ptraj_pred_invalid(vec_invalid==0) = 1-pred_vec_invalid(vec_invalid==0);           
    [r_val, p_val] = corrcoef(ptraj_pred_invalid(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
    [r_inv, p_inv] = corrcoef(ptraj_pred_invalid(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
    predCorr_inval{dataset}(iSubject) = r_val(1,2)-r_inv(1,2);

    ptraj_pred_nochange = zeros(size(pred_vec_nochange));
    ptraj_pred_nochange(vec_nochange==1) = pred_vec_nochange(vec_nochange==1);
    ptraj_pred_nochange(vec_nochange==0) = 1-pred_vec_nochange(vec_nochange==0);           
    [r_val, p_val] = corrcoef(ptraj_pred_nochange(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
    [r_inv, p_inv] = corrcoef(ptraj_pred_nochange(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
    predCorr_nochange{dataset}(iSubject) = r_val(1,2)-r_inv(1,2);

    ptraj_pred_cue = zeros(size(pred_vec_cue));
    ptraj_pred_cue(vec_cue==1) = pred_vec_cue(vec_cue==1);
    ptraj_pred_cue(vec_cue==0) = 1-pred_vec_cue(vec_cue==0);           
    [r_val, p_val] = corrcoef(ptraj_pred_cue(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
    [r_inv, p_inv] = corrcoef(ptraj_pred_cue(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
    predCorr_cue{dataset}(iSubject) = r_val(1,2)-r_inv(1,2);
end


p_sensDelta=[];
p_cc_biasDelta=[];
p_predCorr_val=[];
p_predCorr_inval=[];
p_predCorr_nochange=[];
p_predCorr_cue=[];


for dataset=1:3
    p_sensDelta=[p_sensDelta; sensDelta{dataset}(:)];
    p_cc_biasDelta=[p_cc_biasDelta; cc_biasDelta{dataset}(:)];
    p_predCorr_val=[p_predCorr_val; predCorr_val{dataset}(:)];
    p_predCorr_inval=[p_predCorr_inval; predCorr_inval{dataset}(:)];
    p_predCorr_nochange=[p_predCorr_nochange; predCorr_nochange{dataset}(:)];
    p_predCorr_cue=[p_predCorr_cue; predCorr_cue{dataset}(:)];
end

% sens
[r,~,p] = bendcorr(p_sensDelta(:),p_predCorr_val(:),0);
corr_sensDelta_valid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_sensDelta(:),p_predCorr_inval(:),0);
corr_sensDelta_invalid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_sensDelta(:),p_predCorr_nochange(:),0);
corr_sensDelta_nochange = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_sensDelta(:),p_predCorr_cue(:),0);
corr_sensDelta_cue = struct('r', r, 'p', p);


% bias
[r,~,p] = bendcorr(p_cc_biasDelta(:),p_predCorr_val(:),0);
corr_cc_biasDelta_valid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_cc_biasDelta(:),p_predCorr_inval(:),0);
corr_cc_biasDelta_invalid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_cc_biasDelta(:),p_predCorr_nochange(:),0);
corr_cc_biasDelta_nochange = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_cc_biasDelta(:),p_predCorr_cue(:),0);
corr_cc_biasDelta_cue = struct('r', r, 'p', p);


%% Plotting

xlim_val_d=[-1,3];
xlim_val_cc=[-2.5,1];
ylim_val=[-1,1];

figure(1)
hold on;

subplot(1,4,1); hold on;
scatter(sensDelta{1}(:),predCorr_val{1}(:),'b','o')
scatter(sensDelta{2}(:),predCorr_val{2}(:),'b','s')
scatter(sensDelta{3}(:),predCorr_val{3}(:),'b','^')
title(sess(iSess)+": delta-sens vs delta-corr-valid, "+newline+"r = " + num2str(round(corr_sensDelta_valid.r,3))+", p = " + num2str(round(corr_sensDelta_valid.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("delta sens"); ylabel("delta DBM 'r' ");

subplot(1,4,2); hold on;
scatter(sensDelta{1}(:),predCorr_inval{1}(:),'k','o')
scatter(sensDelta{2}(:),predCorr_inval{2}(:),'k','s')
scatter(sensDelta{3}(:),predCorr_inval{3}(:),'k','^')
title(sess(iSess)+": delta-sens vs delta-corr-invalid, "+newline+"r = " + num2str(round(corr_sensDelta_invalid.r,3))+", p = " + num2str(round(corr_sensDelta_invalid.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("delta sens"); ylabel("delta DBM 'r' ");

subplot(1,4,3); hold on;
scatter(sensDelta{1}(:),predCorr_nochange{1}(:),'r','o')
scatter(sensDelta{2}(:),predCorr_nochange{2}(:),'r','s')
scatter(sensDelta{3}(:),predCorr_nochange{3}(:),'r','^')
title(sess(iSess)+": delta-sens vs delta-corr-nochange, "+newline+"r = " + num2str(round(corr_sensDelta_nochange.r,3))+", p = " + num2str(round(corr_sensDelta_nochange.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("delta sens"); ylabel("delta DBM 'r' ");    

subplot(1,4,4); hold on;
scatter(sensDelta{1}(:),predCorr_cue{1}(:),'g','o')
scatter(sensDelta{2}(:),predCorr_cue{2}(:),'g','s')
scatter(sensDelta{3}(:),predCorr_cue{3}(:),'g','^')
title(sess(iSess)+": delta-sens vs delta-corr-cue, "+newline+"r = " + num2str(round(corr_sensDelta_cue.r,3))+", p = " + num2str(round(corr_sensDelta_cue.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("delta sens"); ylabel("delta DBM 'r' ");


set(gcf,'position',[50,50,1800,400])
%linkaxes([ax1 ax2 ax3 ax4],'xy')
saveas(gcf,plot_savepath+type+"_delta-sens_rt.png")    

figure(2)
hold on;

subplot(1,4,1); hold on;
scatter(cc_biasDelta{1}(:),predCorr_val{1}(:),'b','o')
scatter(cc_biasDelta{2}(:),predCorr_val{2}(:),'b','s')
scatter(cc_biasDelta{3}(:),predCorr_val{3}(:),'b','^')
title(sess(iSess)+": delta-CC bias vs delta-corr-valid, "+newline+"r = " + num2str(round(corr_cc_biasDelta_valid.r,3))+", p = " + num2str(round(corr_cc_biasDelta_valid.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("delta bias"); ylabel("delta DBM 'r' ");

subplot(1,4,2); hold on;
scatter(cc_biasDelta{1}(:),predCorr_inval{1}(:),'k','o')
scatter(cc_biasDelta{2}(:),predCorr_inval{2}(:),'k','s')
scatter(cc_biasDelta{3}(:),predCorr_inval{3}(:),'k','^')
title(sess(iSess)+": delta-CC bias vs delta-corr-invalid, "+newline+"r = " + num2str(round(corr_cc_biasDelta_invalid.r,3))+", p = " + num2str(round(corr_cc_biasDelta_invalid.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("delta bias"); ylabel("delta DBM 'r' ");

subplot(1,4,3); hold on;
scatter(cc_biasDelta{1}(:),predCorr_nochange{1}(:),'r','o')
scatter(cc_biasDelta{2}(:),predCorr_nochange{2}(:),'r','s')
scatter(cc_biasDelta{3}(:),predCorr_nochange{3}(:),'r','^')
title(sess(iSess)+": delta-CC bias vs delta-corr-nochange, "+newline+"r = " + num2str(round(corr_cc_biasDelta_nochange.r,3))+", p = " + num2str(round(corr_cc_biasDelta_nochange.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("delta bias"); ylabel("delta DBM 'r' ");

subplot(1,4,4); hold on;
scatter(cc_biasDelta{1}(:),predCorr_cue{1}(:),'g','o')
scatter(cc_biasDelta{2}(:),predCorr_cue{2}(:),'g','s')
scatter(cc_biasDelta{3}(:),predCorr_cue{3}(:),'g','^')
title(sess(iSess)+": delta-CC bias vs delta-corr-cue, "+newline+"r = " + num2str(round(corr_cc_biasDelta_cue.r,3))+", p = " + num2str(round(corr_cc_biasDelta_cue.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("delta bias"); ylabel("delta DBM 'r' ");

set(gcf,'position',[50,50,1800,400])
%linkaxes([ax1 ax2 ax3 ax4],'xy')
saveas(gcf,plot_savepath+type+"_delta-cc_bias_rt.png")  

