%% Plots the correlation between A and B where
% A = d'/cc/LR bias 
% B = correlation coefficient of predicted prob. (DBM) and RT/accuracy


%% code
clear;
%input type is only "choiceRA"
type="choiceRA";

plot_savepath="results/corr_pred_normDivideByMeanThenLog_";
load("results/compareGridbest_normDivideByMeanThenLog_"+type+".mat") % Load the DBM results file to be plotted

load("data/Ankita_perf_metrics.mat")
for iSide = 1:2
    
    for iSess = 1:3
        
        for iSubject = 1:26
            
            sensVal(iSess,(iSide-1)*26+iSubject) = mean(d_val(iSubject,(iSide-1)*3+iSess,:));
            sensInv(iSess,(iSide-1)*26+iSubject) = mean(d_inv(iSubject,(iSide-1)*3+iSess,:));
            cc_biasVal(iSess,(iSide-1)*26+iSubject) = cc_val(iSubject,(iSide-1)*3+iSess);
            cc_biasInv(iSess,(iSide-1)*26+iSubject) = cc_inv(iSubject,(iSide-1)*3+iSess);
            lr_biasVal(iSess,(iSide-1)*26+iSubject) = get_lr(d_val(iSubject,(iSide-1)*3+iSess,:), c_val(iSubject,(iSide-1)*3+iSess));
            lr_biasInv(iSess,(iSide-1)*26+iSubject) = get_lr(d_inv(iSubject,(iSide-1)*3+iSess,:), c_inv(iSubject,(iSide-1)*3+iSess));
            
            predCorr_choice(iSess,(iSide-1)*26+iSubject) = rvals_choice{iSide,iSess}(iSubject,4);
            
        end
        
    end
    
    
end

for iSess = 1:3
    % sens
   
    [r,~,p] = bendcorr(sensVal(iSess,:),predCorr_choice(iSess,:),0);
    corr_sensVal_choice = struct('r', r, 'p', p);
   
    [r,~,p] = bendcorr(sensInv(iSess,:),predCorr_choice(iSess,:),0);
    corr_sensInv_choice = struct('r', r, 'p', p);
    
    
    % bias
    
    [r,~,p] = bendcorr(cc_biasVal(iSess,:),predCorr_choice(iSess,:),0);
    corr_cc_biasVal_choice = struct('r', r, 'p', p);
    
    [r,~,p] = bendcorr(cc_biasInv(iSess,:),predCorr_choice(iSess,:),0);
    corr_cc_biasInv_choice = struct('r', r, 'p', p);
    
    
    %% Plotting
    
    
    figure(1)
    hold on;
    
    xlim_val_d=[0,3.5];
    xlim_val_cc=[-1,3];
    ylim_val=[-1,1];
    
    subplot(3,1,iSess); hold on; 
    scatter(sensVal(iSess,1:26),predCorr_choice(iSess,1:26),'b')
    scatter(sensVal(iSess,27:52),predCorr_choice(iSess,27:52),'b','filled')
    title(sess(iSess)+": sens-valid vs corr-choiceRA, "+newline+"r = " + num2str(round(corr_sensVal_choice.r,3))+", p = " + num2str(round(corr_sensVal_choice.p,3)));
    xlim(xlim_val_d); xlabel("sens"); ylabel("DBM 'r' ");
    ylim(ylim_val);
    
    set(gcf,'position',[50,50,600,1500])
    saveas(gcf,plot_savepath+type+"_sensVal_rt.png")
    
    figure(2)
    hold on;

    subplot(3,1,iSess); hold on; 
    scatter(sensInv(iSess,1:26),predCorr_choice(iSess,1:26),'b')
    scatter(sensInv(iSess,27:52),predCorr_choice(iSess,27:52),'b','filled')
    title(sess(iSess)+": sens-invalid vs corr-choiceRA, "+newline+"r = " + num2str(round(corr_sensInv_choice.r,3))+", p = " + num2str(round(corr_sensInv_choice.p,3)));
    xlim(xlim_val_d); xlabel("sens"); ylabel("DBM 'r' ");
    ylim(ylim_val);
    
    set(gcf,'position',[50,50,600,1500])
    saveas(gcf,plot_savepath+type+"_sensInv_rt.png")
    

    figure(3)
    hold on;
    
    subplot(3,1,iSess); hold on; 
    scatter(cc_biasVal(iSess,1:26),predCorr_choice(iSess,1:26),'b')
    scatter(cc_biasVal(iSess,27:52),predCorr_choice(iSess,27:52),'b','filled')
    title(sess(iSess)+": CC*bias-valid vs corr-choice, "+newline+"r = " + num2str(round(corr_cc_biasVal_choice.r,3))+", p = " + num2str(round(corr_cc_biasVal_choice.p,3)));
    xlim(xlim_val_cc); xlabel("bias"); ylabel("DBM 'r' ");
    ylim(ylim_val);
    set(gcf,'position',[50,50,600,1500])    
    
    saveas(gcf,plot_savepath+type+"_cc_biasVal_rt.png")
    
    figure(4)
    hold on;
    
    subplot(3,1,iSess); hold on; 
    scatter(cc_biasInv(iSess,1:26),predCorr_choice(iSess,1:26),'b')
    scatter(cc_biasInv(iSess,27:52),predCorr_choice(iSess,27:52),'b','filled')
    title(sess(iSess)+": CC*bias-invalid vs corr-choice, "+newline+"r = " + num2str(round(corr_cc_biasInv_choice.r,3))+", p = " + num2str(round(corr_cc_biasInv_choice.p,3)));
    xlim(xlim_val_cc); xlabel("bias"); ylabel("DBM 'r' ");
    ylim(ylim_val);
    set(gcf,'position',[50,50,600,1500])    
    saveas(gcf,plot_savepath+type+"_cc_biasInv_rt.png")
    
end
