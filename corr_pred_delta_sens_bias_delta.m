% Same as corr_pred_sens_bias, but taking delta-sens (sens(cued) - sens(uncued)) or delta-bias

%% code
clear;
type="respEvent";

plot_savepath="results/corr_pred-delta_";
load("results/compareGridbest_"+type+".mat") % Load the DBM results file to be plotted

load("data/Ankita_perf_metrics.mat")
load("data/tACS_40Hz_woETrej.mat")

for iSide = 1:2
    if iSide==1
        trial_data=leftPPC.Trials_Info;
    else
        trial_data=rightPPC.Trials_Info;
    end
    for iSess = 1:3
        
        for iSubject = 1:26
            
            sensDelta(iSess,(iSide-1)*26+iSubject) = mean(d_val(iSubject,(iSide-1)*3+iSess,:))-mean(d_inv(iSubject,(iSide-1)*3+iSess,:));
            cc_biasDelta(iSess,(iSide-1)*26+iSubject) = cc_val(iSubject,(iSide-1)*3+iSess)-cc_inv(iSubject,(iSide-1)*3+iSess);
                                   
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
            for block=1:5    
                if(type=="respRA1" || var_type=="respRA2")
                    startTrial=2;
                else
                    startTrial=1;
                end
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
            predCorr_val(iSess,(iSide-1)*26+iSubject) = r_val(1,2)-r_inv(1,2);
            
            ptraj_pred_invalid = zeros(size(pred_vec_invalid));
            ptraj_pred_invalid(vec_invalid==1) = pred_vec_invalid(vec_invalid==1);
            ptraj_pred_invalid(vec_invalid==0) = 1-pred_vec_invalid(vec_invalid==0);           
            [r_val, p_val] = corrcoef(ptraj_pred_invalid(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
            [r_inv, p_inv] = corrcoef(ptraj_pred_invalid(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
            predCorr_inval(iSess,(iSide-1)*26+iSubject) = r_val(1,2)-r_inv(1,2);
            
            ptraj_pred_nochange = zeros(size(pred_vec_nochange));
            ptraj_pred_nochange(vec_nochange==1) = pred_vec_nochange(vec_nochange==1);
            ptraj_pred_nochange(vec_nochange==0) = 1-pred_vec_nochange(vec_nochange==0);           
            [r_val, p_val] = corrcoef(ptraj_pred_nochange(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
            [r_inv, p_inv] = corrcoef(ptraj_pred_nochange(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
            predCorr_nochange(iSess,(iSide-1)*26+iSubject) = r_val(1,2)-r_inv(1,2);

            ptraj_pred_cue = zeros(size(pred_vec_cue));
            ptraj_pred_cue(vec_cue==1) = pred_vec_cue(vec_cue==1);
            ptraj_pred_cue(vec_cue==0) = 1-pred_vec_cue(vec_cue==0);           
            [r_val, p_val] = corrcoef(ptraj_pred_cue(rt~=0 & valid_init==1), rt(rt~=0 & valid_init==1));
            [r_inv, p_inv] = corrcoef(ptraj_pred_cue(rt~=0 & invalid_init==1), rt(rt~=0 & invalid_init==1));
            predCorr_cue(iSess,(iSide-1)*26+iSubject) = r_val(1,2)-r_inv(1,2);

            
        end
        
    end
    
    
end

for iSess = 1:3
    % sens
    [r,~,p] = bendcorr(sensDelta(iSess,:),predCorr_val(iSess,:),0);
    corr_sensDelta_valid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(sensDelta(iSess,:),predCorr_inval(iSess,:),0);
    corr_sensDelta_invalid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(sensDelta(iSess,:),predCorr_nochange(iSess,:),0);
    corr_sensDelta_nochange = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(sensDelta(iSess,:),predCorr_cue(iSess,:),0);
    corr_sensDelta_cue = struct('r', r, 'p', p);
   
    
    % bias
    [r,~,p] = bendcorr(cc_biasDelta(iSess,:),predCorr_val(iSess,:),0);
    corr_cc_biasDelta_valid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(cc_biasDelta(iSess,:),predCorr_inval(iSess,:),0);
    corr_cc_biasDelta_invalid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(cc_biasDelta(iSess,:),predCorr_nochange(iSess,:),0);
    corr_cc_biasDelta_nochange = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(cc_biasDelta(iSess,:),predCorr_cue(iSess,:),0);
    corr_cc_biasDelta_cue = struct('r', r, 'p', p);
    
    
    %% Plotting
    
    xlim_val_d=[-1,3];
    xlim_val_cc=[-2,1];
    ylim_val=[-1,1];
    
    figure(1)
    hold on;
    
    subplot(3,4,1+(iSess-1)*4); hold on;
    scatter(sensDelta(iSess,1:26),predCorr_val(iSess,1:26),'b')
    scatter(sensDelta(iSess,27:52),predCorr_val(iSess,27:52),'b','filled')
    title(sess(iSess)+": delta-sens vs delta-corr-valid, "+newline+"r = " + num2str(round(corr_sensDelta_valid.r,3))+", p = " + num2str(round(corr_sensDelta_valid.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("delta sens"); ylabel("delta DBM 'r' ");
    
    subplot(3,4,2+(iSess-1)*4); hold on;
    scatter(sensDelta(iSess,1:26),predCorr_inval(iSess,1:26),'k')
    scatter(sensDelta(iSess,27:52),predCorr_inval(iSess,27:52),'k','filled')
    title(sess(iSess)+": delta-sens vs delta-corr-invalid, "+newline+"r = " + num2str(round(corr_sensDelta_invalid.r,3))+", p = " + num2str(round(corr_sensDelta_invalid.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("delta sens"); ylabel("delta DBM 'r' ");

    subplot(3,4,3+(iSess-1)*4); hold on;
    scatter(sensDelta(iSess,1:26),predCorr_nochange(iSess,1:26),'r')
    scatter(sensDelta(iSess,27:52),predCorr_nochange(iSess,27:52),'r','filled')
    title(sess(iSess)+": delta-sens vs delta-corr-nochange, "+newline+"r = " + num2str(round(corr_sensDelta_nochange.r,3))+", p = " + num2str(round(corr_sensDelta_nochange.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("delta sens"); ylabel("delta DBM 'r' ");    
    
    subplot(3,4,4+(iSess-1)*4); hold on;
    scatter(sensDelta(iSess,1:26),predCorr_cue(iSess,1:26),'g')
    scatter(sensDelta(iSess,27:52),predCorr_cue(iSess,27:52),'g','filled')
    title(sess(iSess)+": delta-sens vs delta-corr-cue, "+newline+"r = " + num2str(round(corr_sensDelta_cue.r,3))+", p = " + num2str(round(corr_sensDelta_cue.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("delta sens"); ylabel("delta DBM 'r' ");
    
   
    set(gcf,'position',[50,50,1800,1000])
    %linkaxes([ax1 ax2 ax3 ax4],'xy')
    saveas(gcf,plot_savepath+type+"_delta-sens_rt.png")    
        
    figure(2)
    hold on;
    
    subplot(3,4,1+(iSess-1)*4); hold on;
    scatter(cc_biasDelta(iSess,1:26),predCorr_val(iSess,1:26),'b')
    scatter(cc_biasDelta(iSess,27:52),predCorr_val(iSess,27:52),'b','filled')
    title(sess(iSess)+": delta-CC bias vs delta-corr-valid, "+newline+"r = " + num2str(round(corr_cc_biasDelta_valid.r,3))+", p = " + num2str(round(corr_cc_biasDelta_valid.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("delta bias"); ylabel("delta DBM 'r' ");
    
    subplot(3,4,2+(iSess-1)*4); hold on;
    scatter(cc_biasDelta(iSess,1:26),predCorr_inval(iSess,1:26),'k')
    scatter(cc_biasDelta(iSess,27:52),predCorr_inval(iSess,27:52),'k','filled')
    title(sess(iSess)+": delta-CC bias vs delta-corr-invalid, "+newline+"r = " + num2str(round(corr_cc_biasDelta_invalid.r,3))+", p = " + num2str(round(corr_cc_biasDelta_invalid.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("delta bias"); ylabel("delta DBM 'r' ");
    
    subplot(3,4,3+(iSess-1)*4); hold on;
    scatter(cc_biasDelta(iSess,1:26),predCorr_nochange(iSess,1:26),'r')
    scatter(cc_biasDelta(iSess,27:52),predCorr_nochange(iSess,27:52),'r','filled')
    title(sess(iSess)+": delta-CC bias vs delta-corr-nochange, "+newline+"r = " + num2str(round(corr_cc_biasDelta_nochange.r,3))+", p = " + num2str(round(corr_cc_biasDelta_nochange.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("delta bias"); ylabel("delta DBM 'r' ");
    
    subplot(3,4,4+(iSess-1)*4); hold on;
    scatter(cc_biasDelta(iSess,1:26),predCorr_cue(iSess,1:26),'g')
    scatter(cc_biasDelta(iSess,27:52),predCorr_cue(iSess,27:52),'g','filled')
    title(sess(iSess)+": delta-CC bias vs delta-corr-cue, "+newline+"r = " + num2str(round(corr_cc_biasDelta_cue.r,3))+", p = " + num2str(round(corr_cc_biasDelta_cue.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("delta bias"); ylabel("delta DBM 'r' ");
    
    set(gcf,'position',[50,50,1800,1000])
    %linkaxes([ax1 ax2 ax3 ax4],'xy')
    saveas(gcf,plot_savepath+type+"_delta-cc_bias_rt.png")  
    
   
   
end
