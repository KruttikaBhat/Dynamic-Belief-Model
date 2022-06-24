clear;
%Event based logic: "respEvent" ; RA1 logic: "respRA1" ; RA2 logic: "respRA2"
type="respEvent";

plot_savepath="results/corr_pred_";
load("results/compareGridbest_"+type+".mat") % Load the DBM results file to be plotted


for iSide = 1:2
    for iSess = 1:3
        for iSubject = 1:26
            predCorr_val(iSess,(iSide-1)*26+iSubject) = rvals_valid{iSide,iSess}(iSubject,4);
            predCorr_inval(iSess,(iSide-1)*26+iSubject) = rvals_invalid{iSide,iSess}(iSubject,4);
        end 
    end
end


for iSess = 1:3    
    %% Plotting
    [r,~,p] = bendcorr(predCorr_val(iSess,:),predCorr_inval(iSess,:),0);
    corr_invVal = struct('r', r, 'p', p);
    
    figure(1)
    hold on;
    
    xlim_val=[-1,1];
    %ylim_val=[-1,0.3];
    ylim_val=[-1,1];
    
    ax1=subplot(3,1,iSess); hold on; 
    scatter(predCorr_val(iSess,1:26),predCorr_inval(iSess,1:26),'b')
    scatter(predCorr_val(iSess,27:52),predCorr_inval(iSess,27:52),'b','filled')
    title(sess(iSess)+": corr-valid vs corr-invalid, "+newline+"r = " + num2str(round(corr_invVal.r,3))+", p = " + num2str(round(corr_invVal.p,3)));
    xlim(xlim_val); 
    ylim(ylim_val);
    xlabel("DBM 'r' valid"); ylabel("DBM 'r' invalid");
    %xticks(0:1:4);
    %yticks(-1:0.5:0.5);
  
    set(gcf,'position',[50,50,600,1500])  
    %linkaxes([ax1 ax2 ax3 ax4],'xy')
    saveas(gcf,plot_savepath+type+"_invVal_DBMr.png")
    

    
end
