% Same as corr_pred_sens_bias, but taking delta-sens (sens(cued) - sens(uncued)) or delta-bias

%% code
clear;
%Event based logic: "respEvent" ; RA1 logic: "respRA1" ; RA2 logic: "respRA2"
type="respEvent";

plot_savepath="results/corr_pred_normDivideByMeanThenLog_";
load("results/compareGridbest_normDivideByMeanThenLog_"+type+".mat") % Load the DBM results file to be plotted

load("data/Ankita_perf_metrics.mat")
for iSide = 1:2
    for iSess = 1:3
        for iSubject = 1:26
            
            sensDelta(iSess,(iSide-1)*26+iSubject) = mean(d_val(iSubject,(iSide-1)*3+iSess,:))-mean(d_inv(iSubject,(iSide-1)*3+iSess,:));
            cc_biasDelta(iSess,(iSide-1)*26+iSubject) = cc_val(iSubject,(iSide-1)*3+iSess)-cc_inv(iSubject,(iSide-1)*3+iSess);
            lr_biasDelta(iSess,(iSide-1)*26+iSubject) = get_lr(d_val(iSubject,(iSide-1)*3+iSess,:),c_val(iSubject,(iSide-1)*3+iSess))-get_lr(d_inv(iSubject,(iSide-1)*3+iSess,:),c_val(iSubject,(iSide-1)*3+iSess));            
          
            predCorr_val(iSess,(iSide-1)*26+iSubject) = rvals_valid{iSide,iSess}(iSubject,4);
            predCorr_inval(iSess,(iSide-1)*26+iSubject) = rvals_invalid{iSide,iSess}(iSubject,4);
            predCorr_nochange(iSess,(iSide-1)*26+iSubject) = rvals_nochange{iSide,iSess}(iSubject,4);
            predCorr_cue(iSess,(iSide-1)*26+iSubject) = rvals_cue{iSide,iSess}(iSubject,4);
            
            
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

    
    % LR Bias
    [r,~,p] = bendcorr(lr_biasDelta(iSess,:),predCorr_val(iSess,:),0);
    corr_lr_biasDelta_valid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(lr_biasDelta(iSess,:),predCorr_inval(iSess,:),0);
    corr_lr_biasDelta_invalid = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(lr_biasDelta(iSess,:),predCorr_nochange(iSess,:),0);
    corr_lr_biasDelta_nochange = struct('r', r, 'p', p);
    [r,~,p] = bendcorr(lr_biasDelta(iSess,:),predCorr_cue(iSess,:),0);
    corr_lr_biasDelta_cue = struct('r', r, 'p', p);
    
    
    %% Plotting
    
    xlim_val_d=[-1,3];
    xlim_val_cc=[-2,1];
    ylim_val=[-1,1];
    
    figure(1)
    hold on;
    
    subplot(3,4,1+(iSess-1)*4); hold on;
    scatter(sensDelta(iSess,1:26),predCorr_val(iSess,1:26),'b')
    scatter(sensDelta(iSess,27:52),predCorr_val(iSess,27:52),'b','filled')
    title(sess(iSess)+": delta-sens vs corr-valid, "+newline+"r = " + num2str(round(corr_sensDelta_valid.r,3))+", p = " + num2str(round(corr_sensDelta_valid.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");
    
    subplot(3,4,2+(iSess-1)*4); hold on;
    scatter(sensDelta(iSess,1:26),predCorr_inval(iSess,1:26),'k')
    scatter(sensDelta(iSess,27:52),predCorr_inval(iSess,27:52),'k','filled')
    title(sess(iSess)+": delta-sens vs corr-invalid, "+newline+"r = " + num2str(round(corr_sensDelta_invalid.r,3))+", p = " + num2str(round(corr_sensDelta_invalid.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");

    subplot(3,4,3+(iSess-1)*4); hold on;
    scatter(sensDelta(iSess,1:26),predCorr_nochange(iSess,1:26),'r')
    scatter(sensDelta(iSess,27:52),predCorr_nochange(iSess,27:52),'r','filled')
    title(sess(iSess)+": delta-sens vs corr-nochange, "+newline+"r = " + num2str(round(corr_sensDelta_nochange.r,3))+", p = " + num2str(round(corr_sensDelta_nochange.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");    
    
    subplot(3,4,4+(iSess-1)*4); hold on;
    scatter(sensDelta(iSess,1:26),predCorr_cue(iSess,1:26),'g')
    scatter(sensDelta(iSess,27:52),predCorr_cue(iSess,27:52),'g','filled')
    title(sess(iSess)+": delta-sens vs corr-cue, "+newline+"r = " + num2str(round(corr_sensDelta_cue.r,3))+", p = " + num2str(round(corr_sensDelta_cue.p,3)));
    xlim(xlim_val_d); 
    ylim(ylim_val);
    xlabel("sens"); ylabel("DBM 'r' ");
    
   
    set(gcf,'position',[50,50,1800,1000])
    %linkaxes([ax1 ax2 ax3 ax4],'xy')
    saveas(gcf,plot_savepath+type+"_delta-sens_rt.png")    
        
    figure(2)
    hold on;
    
    subplot(3,4,1+(iSess-1)*4); hold on;
    scatter(cc_biasDelta(iSess,1:26),predCorr_val(iSess,1:26),'b')
    scatter(cc_biasDelta(iSess,27:52),predCorr_val(iSess,27:52),'b','filled')
    title(sess(iSess)+": delta-CC bias vs corr-valid, "+newline+"r = " + num2str(round(corr_cc_biasDelta_valid.r,3))+", p = " + num2str(round(corr_cc_biasDelta_valid.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
    
    subplot(3,4,2+(iSess-1)*4); hold on;
    scatter(cc_biasDelta(iSess,1:26),predCorr_inval(iSess,1:26),'k')
    scatter(cc_biasDelta(iSess,27:52),predCorr_inval(iSess,27:52),'k','filled')
    title(sess(iSess)+": delta-CC bias vs corr-invalid, "+newline+"r = " + num2str(round(corr_cc_biasDelta_invalid.r,3))+", p = " + num2str(round(corr_cc_biasDelta_invalid.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
    
    subplot(3,4,3+(iSess-1)*4); hold on;
    scatter(cc_biasDelta(iSess,1:26),predCorr_nochange(iSess,1:26),'r')
    scatter(cc_biasDelta(iSess,27:52),predCorr_nochange(iSess,27:52),'r','filled')
    title(sess(iSess)+": delta-CC bias vs corr-nochange, "+newline+"r = " + num2str(round(corr_cc_biasDelta_nochange.r,3))+", p = " + num2str(round(corr_cc_biasDelta_nochange.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
    
    subplot(3,4,4+(iSess-1)*4); hold on;
    scatter(cc_biasDelta(iSess,1:26),predCorr_cue(iSess,1:26),'g')
    scatter(cc_biasDelta(iSess,27:52),predCorr_cue(iSess,27:52),'g','filled')
    title(sess(iSess)+": delta-CC bias vs corr-cue, "+newline+"r = " + num2str(round(corr_cc_biasDelta_cue.r,3))+", p = " + num2str(round(corr_cc_biasDelta_cue.p,3)));
    xlim(xlim_val_cc); 
    ylim(ylim_val);
    xlabel("bias"); ylabel("DBM 'r' ");
    
    set(gcf,'position',[50,50,1800,1000])
    %linkaxes([ax1 ax2 ax3 ax4],'xy')
    saveas(gcf,plot_savepath+type+"_delta-cc_bias_rt.png")  
    
   
   
end
