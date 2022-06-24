%% Plots the correlation between A and B where
% A = d'/cc/LR bias 
% B = correlation coefficient of predicted prob. (DBM) and RT/accuracy


%% code
clear;
type="respEvent";

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

            predCorr_val(iSess,(iSide-1)*26+iSubject) = rvals_valid{iSide,iSess}(iSubject,4);
            predCorr_inval(iSess,(iSide-1)*26+iSubject) = rvals_invalid{iSide,iSess}(iSubject,4);
            predCorr_nochange(iSess,(iSide-1)*26+iSubject) = rvals_nochange{iSide,iSess}(iSubject,4);
            predCorr_cue(iSess,(iSide-1)*26+iSubject) = rvals_cue{iSide,iSess}(iSubject,4);
             
            
        end
        
    end
    
    
end

for iSess = 1:3
    % sens
    [r,~,p] = bendcorr(sensVal(iSess,:),predCorr_val(iSess,:),0);
    corr_sensVal_valid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(sensVal(iSess,:),predCorr_inval(iSess,:),0);
    corr_sensVal_invalid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(sensVal(iSess,:),predCorr_nochange(iSess,:),0);
    corr_sensVal_nochange = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(sensVal(iSess,:),predCorr_cue(iSess,:),0);
    corr_sensVal_cue = struct('r', r, 'p', p);
    
    [r,~,p] = bendcorr(sensInv(iSess,:),predCorr_val(iSess,:),0);
    corr_sensInv_valid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(sensInv(iSess,:),predCorr_inval(iSess,:),0);
    corr_sensInv_invalid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(sensInv(iSess,:),predCorr_nochange(iSess,:),0);
    corr_sensInv_nochange = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(sensInv(iSess,:),predCorr_cue(iSess,:),0);
    corr_sensInv_cue = struct('r', r, 'p', p);
    
    
    % bias
    [r,~,p] = bendcorr(cc_biasVal(iSess,:),predCorr_val(iSess,:),0);
    corr_cc_biasVal_valid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(cc_biasVal(iSess,:),predCorr_inval(iSess,:),0);
    corr_cc_biasVal_invalid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(cc_biasVal(iSess,:),predCorr_nochange(iSess,:),0);
    corr_cc_biasVal_nochange = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(cc_biasVal(iSess,:),predCorr_cue(iSess,:),0);
    corr_cc_biasVal_cue = struct('r', r, 'p', p);
    
    [r,~,p] = bendcorr(cc_biasInv(iSess,:),predCorr_val(iSess,:),0);
    corr_cc_biasInv_valid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(cc_biasInv(iSess,:),predCorr_inval(iSess,:),0);
    corr_cc_biasInv_invalid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(cc_biasInv(iSess,:),predCorr_nochange(iSess,:),0);
    corr_cc_biasInv_nochange = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(cc_biasInv(iSess,:),predCorr_cue(iSess,:),0);
    corr_cc_biasInv_cue = struct('r', r, 'p', p);
    
    
    % LR Bias
    [r,~,p] = bendcorr(lr_biasVal(iSess,:),predCorr_val(iSess,:),0);
    corr_lr_biasVal_valid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(lr_biasVal(iSess,:),predCorr_inval(iSess,:),0);
    corr_lr_biasVal_invalid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(lr_biasVal(iSess,:),predCorr_nochange(iSess,:),0);
    corr_lr_biasVal_nochange = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(lr_biasVal(iSess,:),predCorr_cue(iSess,:),0);
    corr_lr_biasVal_cue = struct('r', r, 'p', p);
    
    [r,~,p] = bendcorr(lr_biasInv(iSess,:),predCorr_val(iSess,:),0);
    corr_lr_biasInv_valid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(lr_biasInv(iSess,:),predCorr_inval(iSess,:),0);
    corr_lr_biasInv_invalid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(lr_biasInv(iSess,:),predCorr_nochange(iSess,:),0);
    corr_lr_biasInv_nochange = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(lr_biasInv(iSess,:),predCorr_cue(iSess,:),0);
    corr_lr_biasInv_cue = struct('r', r, 'p', p);
  
    
    %% Plotting
    
    
    figure(1)
    hold on;
    
    xlim_val_d=[0,3.5];
    xlim_val_cc=[-1,3];
    %ylim_val=[-1,0.3];
    ylim_val=[-1,1];
    
    ax1=subplot(3,4,1+(iSess-1)*4); hold on; 
    scatter(sensVal(iSess,1:26),predCorr_val(iSess,1:26),'b')
    scatter(sensVal(iSess,27:52),predCorr_val(iSess,27:52),'b','filled')
    title(sess(iSess)+": sens-valid vs corr-valid, "+newline+"r = " + num2str(round(corr_sensVal_valid.r,3))+", p = " + num2str(round(corr_sensVal_valid.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");
    %xticks(0:1:4);
    %yticks(-1:0.5:0.5);
    
    ax2=subplot(3,4,2+(iSess-1)*4); hold on; 
    scatter(sensVal(iSess,1:26),predCorr_inval(iSess,1:26),'k')
    scatter(sensVal(iSess,27:52),predCorr_inval(iSess,27:52),'k','filled')
    title(sess(iSess)+": sens-valid vs corr-invalid, "+newline+"r = " + num2str(round(corr_sensVal_invalid.r,3))+", p = " + num2str(round(corr_sensVal_invalid.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");
    %xticks(0:1:4);
    %yticks(-1:0.5:0.5);
    
    ax3=subplot(3,4,3+(iSess-1)*4); hold on; 
    scatter(sensVal(iSess,1:26),predCorr_nochange(iSess,1:26),'r')
    scatter(sensVal(iSess,27:52),predCorr_nochange(iSess,27:52),'r','filled')
    title(sess(iSess)+": sens-valid vs corr-nochange, "+newline+"r = " + num2str(round(corr_sensVal_nochange.r,3))+", p = " + num2str(round(corr_sensVal_nochange.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");
    %xticks(0:1:4);
    %yticks(-1:0.5:0.5);
    
    ax4=subplot(3,4,4+(iSess-1)*4); hold on; 
    scatter(sensVal(iSess,1:26),predCorr_cue(iSess,1:26),'g')
    scatter(sensVal(iSess,27:52),predCorr_cue(iSess,27:52),'g','filled')
    title(sess(iSess)+": sens-valid vs corr-cue, "+newline+"r = " + num2str(round(corr_sensVal_cue.r,3))+", p = " + num2str(round(corr_sensVal_cue.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");
    %xticks(0:1:4);
    %yticks(-1:0.5:0.5);
   
    set(gcf,'position',[50,50,1800,1000])
    %linkaxes([ax1 ax2 ax3 ax4],'xy')
    saveas(gcf,plot_savepath+type+"_sensVal_rt.png")
    
    figure(2)
    hold on;
    
    
    ax1=subplot(3,4,1+(iSess-1)*4); hold on; 
    scatter(sensInv(iSess,1:26),predCorr_val(iSess,1:26),'b')
    scatter(sensInv(iSess,27:52),predCorr_val(iSess,27:52),'b','filled')
    title(sess(iSess)+": sens-invalid vs corr-valid, "+newline+"r = " + num2str(round(corr_sensInv_valid.r,3))+", p = " + num2str(round(corr_sensInv_valid.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);
    
    ax2=subplot(3,4,2+(iSess-1)*4); hold on; 
    scatter(sensInv(iSess,1:26),predCorr_inval(iSess,1:26),'k')
    scatter(sensInv(iSess,27:52),predCorr_inval(iSess,27:52),'k','filled')
    title(sess(iSess)+": sens-invalid vs corr-invalid, "+newline+"r = " + num2str(round(corr_sensInv_invalid.r,3))+", p = " + num2str(round(corr_sensInv_invalid.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);
    
    ax3=subplot(3,4,3+(iSess-1)*4); hold on; 
    scatter(sensInv(iSess,1:26),predCorr_nochange(iSess,1:26),'r')
    scatter(sensInv(iSess,27:52),predCorr_nochange(iSess,27:52),'r','filled')
    title(sess(iSess)+": sens-invalid vs corr-nochange, "+newline+"r = " + num2str(round(corr_sensInv_nochange.r,3))+", p = " + num2str(round(corr_sensInv_nochange.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);
    
    ax4=subplot(3,4,4+(iSess-1)*4); hold on; 
    scatter(sensInv(iSess,1:26),predCorr_cue(iSess,1:26),'g')
    scatter(sensInv(iSess,27:52),predCorr_cue(iSess,27:52),'g','filled')
    title(sess(iSess)+": sens-invalid vs corr-cue, "+newline+"r = " + num2str(round(corr_sensInv_cue.r,3))+", p = " + num2str(round(corr_sensInv_cue.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);
    set(gcf,'position',[50,50,1800,1500])
    %linkaxes([ax1 ax2 ax3 ax4],'xy')
    saveas(gcf,plot_savepath+type+"_sensInv_rt.png")
    
    figure(3)
    hold on;
    
    ax1=subplot(3,4,1+(iSess-1)*4); hold on; 
    scatter(cc_biasVal(iSess,1:26),predCorr_val(iSess,1:26),'b')
    scatter(cc_biasVal(iSess,27:52),predCorr_val(iSess,27:52),'b','filled')
    title(sess(iSess)+": CC*bias-valid vs corr-valid, "+newline+"r = " + num2str(round(corr_cc_biasVal_valid.r,3))+", p = " + num2str(round(corr_cc_biasVal_valid.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);

    ax2=subplot(3,4,2+(iSess-1)*4); hold on; 
    scatter(cc_biasVal(iSess,1:26),predCorr_inval(iSess,1:26),'k')
    scatter(cc_biasVal(iSess,27:52),predCorr_inval(iSess,27:52),'k','filled')
    title(sess(iSess)+": CC*bias valid vs corr-invalid, "+newline+"r = " + num2str(round(corr_cc_biasVal_invalid.r,3))+", p = " + num2str(round(corr_cc_biasVal_invalid.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);
    
    ax3=subplot(3,4,3+(iSess-1)*4); hold on; 
    scatter(cc_biasVal(iSess,1:26),predCorr_nochange(iSess,1:26),'r')
    scatter(cc_biasVal(iSess,27:52),predCorr_nochange(iSess,27:52),'r','filled')
    title(sess(iSess)+": CC*bias valid vs corr-nochange, "+newline+"r = " + num2str(round(corr_cc_biasVal_nochange.r,3))+", p = " + num2str(round(corr_cc_biasVal_nochange.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);
    
    ax4=subplot(3,4,4+(iSess-1)*4); hold on; 
    scatter(cc_biasVal(iSess,1:26),predCorr_cue(iSess,1:26),'g')
    scatter(cc_biasVal(iSess,27:52),predCorr_cue(iSess,27:52),'g','filled')
    title(sess(iSess)+": CC*bias valid vs corr-cue, "+newline+"r = " + num2str(round(corr_cc_biasVal_cue.r,3))+", p = " + num2str(round(corr_cc_biasVal_cue.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);
    set(gcf,'position',[50,50,1800,1500])    
    %linkaxes([ax1 ax2 ax3 ax4],'xy')
    saveas(gcf,plot_savepath+type+"_cc_biasVal_rt.png")
    
    figure(4)
    hold on;
    
    ax1=subplot(3,4,1+(iSess-1)*4); hold on; 
    scatter(cc_biasInv(iSess,1:26),predCorr_val(iSess,1:26),'b')
    scatter(cc_biasInv(iSess,27:52),predCorr_val(iSess,27:52),'b','filled')
    title(sess(iSess)+": CC*bias-invalid vs corr-valid, "+newline+"r = " + num2str(round(corr_cc_biasInv_valid.r,3))+", p = " + num2str(round(corr_cc_biasInv_valid.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);

    ax2=subplot(3,4,2+(iSess-1)*4); hold on; 
    scatter(cc_biasInv(iSess,1:26),predCorr_inval(iSess,1:26),'k')
    scatter(cc_biasInv(iSess,27:52),predCorr_inval(iSess,27:52),'k','filled')
    title(sess(iSess)+": CC*bias invalid vs corr-invalid, "+newline+"r = " + num2str(round(corr_cc_biasInv_invalid.r,3))+", p = " + num2str(round(corr_cc_biasInv_invalid.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);
    
    ax3=subplot(3,4,3+(iSess-1)*4); hold on; 
    scatter(cc_biasInv(iSess,1:26),predCorr_nochange(iSess,1:26),'r')
    scatter(cc_biasInv(iSess,27:52),predCorr_nochange(iSess,27:52),'r','filled')
    title(sess(iSess)+": CC*bias invalid vs corr-nochange, "+newline+"r = " + num2str(round(corr_cc_biasInv_nochange.r,3))+", p = " + num2str(round(corr_cc_biasInv_nochange.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);
    
    ax4=subplot(3,4,4+(iSess-1)*4); hold on; 
    scatter(cc_biasInv(iSess,1:26),predCorr_cue(iSess,1:26),'g')
    scatter(cc_biasInv(iSess,27:52),predCorr_cue(iSess,27:52),'g','filled')
    title(sess(iSess)+": CC*bias invalid vs corr-cue, "+newline+"r = " + num2str(round(corr_cc_biasInv_cue.r,3))+", p = " + num2str(round(corr_cc_biasInv_cue.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
%     xticks(0:1:4);
%     yticks(-1:0.5:0.5);
    set(gcf,'position',[50,50,1800,1500])  
    %linkaxes([ax1 ax2 ax3 ax4],'xy')
    saveas(gcf,plot_savepath+type+"_cc_biasInv_rt.png")
    
%     figure(5)
%     hold on;
%     subplot(3,3,1+(iSess-1)*3); hold on;
%     scatter(lr_biasVal(iSess,1:26),predCorr_isval(iSess,1:26),'b')
%     scatter(lr_biasVal(iSess,27:52),predCorr_isval(iSess,27:52),'b','filled')
%     
%     title(sess(iSess)+": lr bias-valid vs corr-isvalid, r = " + num2str(round(corr_lr_biasVal_isvalid.r,3))+", p = " + num2str(round(corr_lr_biasVal_isvalid.p,3)));
%     xlim([-0.15,3.3]); xlabel("bias"); ylabel("DBM 'r' ");
%     ylim([-0.2,0.8])
%     subplot(3,3,2+(iSess-1)*3); hold on;
%     scatter(lr_biasVal(iSess,1:26),predCorr_isleft(iSess,1:26),'k')
%     scatter(lr_biasVal(iSess,27:52),predCorr_isleft(iSess,27:52),'k','filled')
%     
%     title(sess(iSess)+": lr bias-valid vs corr-left, r = " + num2str(round(corr_lr_biasVal_isleft.r,3))+", p = " + num2str(round(corr_lr_biasVal_isleft.p,3)));
%     xlim([-0.15,3.3]); xlabel("bias"); ylabel("DBM 'r' ");
%     ylim([-0.2,0.8])
%     subplot(3,3,3+(iSess-1)*3); hold on;
%     scatter(lr_biasVal(iSess,1:26),predCorr_ischange(iSess,1:26),'r')
%     scatter(lr_biasVal(iSess,27:52),predCorr_ischange(iSess,27:52),'r','filled')
%     
%     title(sess(iSess)+": lr bias-valid vs corr-ischange, r = " + num2str(round(corr_lr_biasVal_ischange.r,3))+", p = " + num2str(round(corr_lr_biasVal_ischange.p,3)));
%     xlim([-0.15,3.3]); xlabel("bias"); ylabel("DBM 'r' ");
%     ylim([-0.2,0.8])
%     
%     set(gcf,'position',[50,400,1600,1500])
%     
%     figure(6)
%     hold on;
%     subplot(3,3,1+(iSess-1)*3); hold on;
%     scatter(lr_biasInv(iSess,1:26),predCorr_isval(iSess,1:26),'b')
%     
%     scatter(lr_biasInv(iSess,27:52),predCorr_isval(iSess,27:52),'b','filled')
%     title(sess(iSess)+": lr bias-invalid vs corr-isvalid, r = " + num2str(round(corr_lr_biasInv_isvalid.r,3))+", p = " + num2str(round(corr_lr_biasInv_isvalid.p,3)));
%     xlim([-0.15,3.3]); xlabel("bias"); ylabel("DBM 'r' ");
%     ylim([-0.2,0.8])
%     subplot(3,3,2+(iSess-1)*3); hold on;
%     scatter(lr_biasInv(iSess,1:26),predCorr_isleft(iSess,1:26),'k')
%     scatter(lr_biasInv(iSess,27:52),predCorr_isleft(iSess,27:52),'k','filled')
%     
%     title(sess(iSess)+": lr bias-invalid vs corr-left, r = " + num2str(round(corr_lr_biasInv_isleft.r,3))+", p = " + num2str(round(corr_lr_biasInv_isleft.p,3)));
%     xlim([-0.15,3.3]); xlabel("bias"); ylabel("DBM 'r' ");
%     ylim([-0.2,0.8])
%     subplot(3,3,3+(iSess-1)*3); hold on;
%     scatter(lr_biasInv(iSess,1:26),predCorr_ischange(iSess,1:26),'r')
%     scatter(lr_biasInv(iSess,27:52),predCorr_ischange(iSess,27:52),'r','filled')
%     
%     title(sess(iSess)+": lr bias-invalid vs corr-ischange, r = " + num2str(round(corr_lr_biasInv_ischange.r,3))+", p = " + num2str(round(corr_lr_biasInv_ischange.p,3)));
%     xlim([-0.15,3.3]); xlabel("bias"); ylabel("DBM 'r' ");
%     ylim([-0.2,0.8])
%     
%     set(gcf,'position',[10,10,1600,1500])
    
end
