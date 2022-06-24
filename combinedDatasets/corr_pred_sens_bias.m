%% Plots the correlation between A and B where
% A = d'/cc/LR bias 
% B = correlation coefficient of predicted prob. (DBM) and RT/accuracy


%% code
clear;
type="respRA2";
sess=["Sham"];
plot_savepath="results/corr_pred_";

%get Ankita's data
load("../DBM_fit_Ankita_data/results/compareGridbest_"+type+".mat") % Load the DBM results file to be plotted
load("../DBM_fit_Ankita_data/data/Ankita_perf_metrics.mat")
iSess=1;
dataset=1;

for iSide = 1:2
    for iSubject = 1:26            
        sensVal{dataset}((iSide-1)*26+iSubject) = mean(d_val(iSubject,(iSide-1)*3+iSess,:));
        sensInv{dataset}((iSide-1)*26+iSubject) = mean(d_inv(iSubject,(iSide-1)*3+iSess,:));
        cc_biasVal{dataset}((iSide-1)*26+iSubject) = cc_val(iSubject,(iSide-1)*3+iSess);
        cc_biasInv{dataset}((iSide-1)*26+iSubject) = cc_inv(iSubject,(iSide-1)*3+iSess);
        
        predCorr_val{dataset}((iSide-1)*26+iSubject) = rvals_valid{iSide,iSess}(iSubject,4);
        predCorr_inval{dataset}((iSide-1)*26+iSubject) = rvals_invalid{iSide,iSess}(iSubject,4);
        predCorr_nochange{dataset}((iSide-1)*26+iSubject) = rvals_nochange{iSide,iSess}(iSubject,4);
        predCorr_cue{dataset}((iSide-1)*26+iSubject) = rvals_cue{iSide,iSess}(iSubject,4);

    end
end
  
%get Sanjna's data
load("../DBM_fit_Sanjna_data/results/compareGridbest_"+type+".mat") % Load the DBM results file to be plotted
load("../DBM_fit_Sanjna_data/data/Sanjna_perf_metrics.mat")

iSess=1;
dataset=2;

for iSubject = 1:28
    sensVal{dataset}(iSubject) = mean(d_val(iSubject,iSess,:));
    sensInv{dataset}(iSubject) = mean(d_inv(iSubject,iSess,:));
    cc_biasVal{dataset}(iSubject) = cc_val(iSubject,iSess);
    cc_biasInv{dataset}(iSubject) = cc_inv(iSubject,iSess);
    
    predCorr_val{dataset}(iSubject) = rvals_valid{iSess}(iSubject,4);
    predCorr_inval{dataset}(iSubject) = rvals_invalid{iSess}(iSubject,4);
    predCorr_nochange{dataset}(iSubject) = rvals_nochange{iSess}(iSubject,4);
    predCorr_cue{dataset}(iSubject) = rvals_cue{iSess}(iSubject,4);
end


%get Varsha's data
load("../DBM_fit_Varsha_data/results/10Blocks/compareGridbest_"+type+".mat") % Load the DBM results file to be plotted
load("../DBM_fit_Varsha_data/data/Varsha_10Blocks_perf_metrics.mat")

iSess=1;
dataset=3;

for iSubject = 1:22
    sensVal{dataset}(iSubject) = mean(d_val(iSubject,iSess,:));
    sensInv{dataset}(iSubject) = mean(d_inv(iSubject,iSess,:));
    cc_biasVal{dataset}(iSubject) = cc_val(iSubject,iSess);
    cc_biasInv{dataset}(iSubject) = cc_inv(iSubject,iSess);
    
    predCorr_val{dataset}(iSubject) = rvals_valid{iSess}(iSubject,4);
    predCorr_inval{dataset}(iSubject) = rvals_invalid{iSess}(iSubject,4);
    predCorr_nochange{dataset}(iSubject) = rvals_nochange{iSess}(iSubject,4);
    predCorr_cue{dataset}(iSubject) = rvals_cue{iSess}(iSubject,4);
end

p_sensVal=[];
p_cc_biasVal=[];
p_sensInv=[];
p_cc_biasInv=[];
p_predCorr_val=[];
p_predCorr_inval=[];
p_predCorr_nochange=[];
p_predCorr_cue=[];


for dataset=1:3
    p_sensVal=[p_sensVal; sensVal{dataset}(:)];
    p_cc_biasVal=[p_cc_biasVal; cc_biasVal{dataset}(:)];
    p_sensInv=[p_sensInv; sensInv{dataset}(:)];
    p_cc_biasInv=[p_cc_biasInv; cc_biasInv{dataset}(:)];
    p_predCorr_val=[p_predCorr_val; predCorr_val{dataset}(:)];
    p_predCorr_inval=[p_predCorr_inval; predCorr_inval{dataset}(:)];
    p_predCorr_nochange=[p_predCorr_nochange; predCorr_nochange{dataset}(:)];
    p_predCorr_cue=[p_predCorr_cue; predCorr_cue{dataset}(:)];
end

% sens
[r,~,p] = bendcorr(p_sensVal(:),p_predCorr_val(:),0);
corr_sensVal_valid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_sensVal(:),p_predCorr_inval(:),0);
corr_sensVal_invalid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_sensVal(:),p_predCorr_nochange(:),0);
corr_sensVal_nochange = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_sensVal(:),p_predCorr_cue(:),0);
corr_sensVal_cue = struct('r', r, 'p', p);

[r,~,p] = bendcorr(p_sensInv(:),p_predCorr_val(:),0);
corr_sensInv_valid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_sensInv(:),p_predCorr_inval(:),0);
corr_sensInv_invalid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_sensInv(:),p_predCorr_nochange(:),0);
corr_sensInv_nochange = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_sensInv(:),p_predCorr_cue(:),0);
corr_sensInv_cue = struct('r', r, 'p', p);


% bias
[r,~,p] = bendcorr(p_cc_biasVal(:),p_predCorr_val(:),0);
corr_cc_biasVal_valid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_cc_biasVal(:),p_predCorr_inval(:),0);
corr_cc_biasVal_invalid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_cc_biasVal(:),p_predCorr_nochange(:),0);
corr_cc_biasVal_nochange = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_cc_biasVal(:),p_predCorr_cue(:),0);
corr_cc_biasVal_cue = struct('r', r, 'p', p);

[r,~,p] = bendcorr(p_cc_biasInv(:),p_predCorr_val(:),0);
corr_cc_biasInv_valid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_cc_biasInv(:),p_predCorr_inval(:),0);
corr_cc_biasInv_invalid = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_cc_biasInv(:),p_predCorr_nochange(:),0);
corr_cc_biasInv_nochange = struct('r', r, 'p', p);
[r,~,p] = bendcorr(p_cc_biasInv(:),p_predCorr_cue(:),0);
corr_cc_biasInv_cue = struct('r', r, 'p', p);


%% Plotting

xlim_val_d=[0,3.5];
xlim_val_cc=[-1,1.5];
ylim_val=[-1,1];

figure(1)
hold on;

subplot(1,4,1); hold on;
scatter(sensVal{1}(:),predCorr_val{1}(:),'b','o')
scatter(sensVal{2}(:),predCorr_val{2}(:),'b','s')
scatter(sensVal{3}(:),predCorr_val{3}(:),'b','^')
title(sess(iSess)+": sens-valid vs corr-valid, "+newline+"r = " + num2str(round(corr_sensVal_valid.r,3))+", p = " + num2str(round(corr_sensVal_valid.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("sens"); ylabel("DBM 'r' ");

subplot(1,4,2); hold on;
scatter(sensVal{1}(:),predCorr_inval{1}(:),'k','o')
scatter(sensVal{2}(:),predCorr_inval{2}(:),'k','s')
scatter(sensVal{3}(:),predCorr_inval{3}(:),'k','^')
title(sess(iSess)+": sens-valid vs corr-invalid, "+newline+"r = " + num2str(round(corr_sensVal_invalid.r,3))+", p = " + num2str(round(corr_sensVal_invalid.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("sens"); ylabel("DBM 'r' ");

subplot(1,4,3); hold on;
scatter(sensVal{1}(:),predCorr_nochange{1}(:),'r','o')
scatter(sensVal{2}(:),predCorr_nochange{2}(:),'r','s')
scatter(sensVal{3}(:),predCorr_nochange{3}(:),'r','^')
title(sess(iSess)+": sens-valid vs corr-nochange, "+newline+"r = " + num2str(round(corr_sensVal_nochange.r,3))+", p = " + num2str(round(corr_sensVal_nochange.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("sens"); ylabel("DBM 'r' ");

subplot(1,4,4); hold on;
scatter(sensVal{1}(:),predCorr_cue{1}(:),'g','o')
scatter(sensVal{2}(:),predCorr_cue{2}(:),'g','s')
scatter(sensVal{3}(:),predCorr_cue{3}(:),'g','^')
title(sess(iSess)+": sens-valid vs corr-cue, "+newline+"r = " + num2str(round(corr_sensVal_cue.r,3))+", p = " + num2str(round(corr_sensVal_cue.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("sens"); ylabel("DBM 'r' ");

set(gcf,'position',[50,50,1800,400])
%linkaxes([ax1 ax2 ax3 ax4],'xy')
saveas(gcf,plot_savepath+type+"_sensVal_rt.png")    


figure(2)
hold on;

subplot(1,4,1); hold on;
scatter(sensInv{1}(:),predCorr_val{1}(:),'b','o')
scatter(sensInv{2}(:),predCorr_val{2}(:),'b','s')
scatter(sensInv{3}(:),predCorr_val{3}(:),'b','^')
title(sess(iSess)+": sens-invalid vs corr-valid, "+newline+"r = " + num2str(round(corr_sensInv_valid.r,3))+", p = " + num2str(round(corr_sensInv_valid.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("sens"); ylabel("DBM 'r' ");

subplot(1,4,2); hold on;
scatter(sensInv{1}(:),predCorr_inval{1}(:),'k','o')
scatter(sensInv{2}(:),predCorr_inval{2}(:),'k','s')
scatter(sensInv{3}(:),predCorr_inval{3}(:),'k','^')
title(sess(iSess)+": sens-invalid vs corr-invalid, "+newline+"r = " + num2str(round(corr_sensInv_invalid.r,3))+", p = " + num2str(round(corr_sensInv_invalid.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("sens"); ylabel("DBM 'r' ");

subplot(1,4,3); hold on;
scatter(sensInv{1}(:),predCorr_nochange{1}(:),'r','o')
scatter(sensInv{2}(:),predCorr_nochange{2}(:),'r','s')
scatter(sensInv{3}(:),predCorr_nochange{3}(:),'r','^')
title(sess(iSess)+": sens-invalid vs corr-nochange, "+newline+"r = " + num2str(round(corr_sensInv_nochange.r,3))+", p = " + num2str(round(corr_sensInv_nochange.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("sens"); ylabel("DBM 'r' ");

subplot(1,4,4); hold on;
scatter(sensInv{1}(:),predCorr_cue{1}(:),'g','o')
scatter(sensInv{2}(:),predCorr_cue{2}(:),'g','s')
scatter(sensInv{3}(:),predCorr_cue{3}(:),'g','^')
title(sess(iSess)+": sens-invalid vs corr-cue, "+newline+"r = " + num2str(round(corr_sensInv_cue.r,3))+", p = " + num2str(round(corr_sensInv_cue.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("sens"); ylabel("DBM 'r' ");

set(gcf,'position',[50,50,1800,400])
%linkaxes([ax1 ax2 ax3 ax4],'xy')
saveas(gcf,plot_savepath+type+"_sensInv_rt.png")   


figure(3)
hold on;

subplot(1,4,1); hold on;
scatter(cc_biasVal{1}(:),predCorr_val{1}(:),'b','o')
scatter(cc_biasVal{2}(:),predCorr_val{2}(:),'b','s')
scatter(cc_biasVal{3}(:),predCorr_val{3}(:),'b','^')
title(sess(iSess)+": CC*bias-valid vs corr-valid, "+newline+"r = " + num2str(round(corr_cc_biasVal_valid.r,3))+", p = " + num2str(round(corr_cc_biasVal_valid.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("bias"); ylabel("DBM 'r' ");

subplot(1,4,2); hold on;
scatter(cc_biasVal{1}(:),predCorr_inval{1}(:),'k','o')
scatter(cc_biasVal{2}(:),predCorr_inval{2}(:),'k','s')
scatter(cc_biasVal{3}(:),predCorr_inval{3}(:),'k','^')
title(sess(iSess)+": CC*bias-valid vs corr-invalid, "+newline+"r = " + num2str(round(corr_cc_biasVal_invalid.r,3))+", p = " + num2str(round(corr_cc_biasVal_invalid.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("bias"); ylabel("DBM 'r' ");

subplot(1,4,3); hold on;
scatter(cc_biasVal{1}(:),predCorr_nochange{1}(:),'r','o')
scatter(cc_biasVal{2}(:),predCorr_nochange{2}(:),'r','s')
scatter(cc_biasVal{3}(:),predCorr_nochange{3}(:),'r','^')
title(sess(iSess)+": CC*bias-valid vs corr-nochange, "+newline+"r = " + num2str(round(corr_cc_biasVal_nochange.r,3))+", p = " + num2str(round(corr_cc_biasVal_nochange.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("bias"); ylabel("DBM 'r' ");

subplot(1,4,4); hold on;
scatter(cc_biasVal{1}(:),predCorr_cue{1}(:),'g','o')
scatter(cc_biasVal{2}(:),predCorr_cue{2}(:),'g','s')
scatter(cc_biasVal{3}(:),predCorr_cue{3}(:),'g','^')
title(sess(iSess)+": CC*bias-valid vs corr-cue, "+newline+"r = " + num2str(round(corr_cc_biasVal_cue.r,3))+", p = " + num2str(round(corr_cc_biasVal_cue.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("bias"); ylabel("DBM 'r' ");

set(gcf,'position',[50,50,1800,400])
%linkaxes([ax1 ax2 ax3 ax4],'xy')
saveas(gcf,plot_savepath+type+"_cc_biasVal_rt.png")

    
figure(4)
hold on;

subplot(1,4,1); hold on;
scatter(cc_biasInv{1}(:),predCorr_val{1}(:),'b','o')
scatter(cc_biasInv{2}(:),predCorr_val{2}(:),'b','s')
scatter(cc_biasInv{3}(:),predCorr_val{3}(:),'b','^')
title(sess(iSess)+": CC*bias-invalid vs corr-valid, "+newline+"r = " + num2str(round(corr_cc_biasInv_valid.r,3))+", p = " + num2str(round(corr_cc_biasInv_valid.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("bias"); ylabel("DBM 'r' ");

subplot(1,4,2); hold on;
scatter(cc_biasInv{1}(:),predCorr_inval{1}(:),'k','o')
scatter(cc_biasInv{2}(:),predCorr_inval{2}(:),'k','s')
scatter(cc_biasInv{3}(:),predCorr_inval{3}(:),'k','^')
title(sess(iSess)+": CC*bias-invalid vs corr-invalid, "+newline+"r = " + num2str(round(corr_cc_biasInv_invalid.r,3))+", p = " + num2str(round(corr_cc_biasInv_invalid.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("bias"); ylabel("DBM 'r' ");

subplot(1,4,3); hold on;
scatter(cc_biasInv{1}(:),predCorr_nochange{1}(:),'r','o')
scatter(cc_biasInv{2}(:),predCorr_nochange{2}(:),'r','s')
scatter(cc_biasInv{3}(:),predCorr_nochange{3}(:),'r','^')
title(sess(iSess)+": CC*bias-invalid vs corr-nochange, "+newline+"r = " + num2str(round(corr_cc_biasInv_nochange.r,3))+", p = " + num2str(round(corr_cc_biasInv_nochange.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("bias"); ylabel("DBM 'r' ");

subplot(1,4,4); hold on;
scatter(cc_biasInv{1}(:),predCorr_cue{1}(:),'g','o')
scatter(cc_biasInv{2}(:),predCorr_cue{2}(:),'g','s')
scatter(cc_biasInv{3}(:),predCorr_cue{3}(:),'g','^')
title(sess(iSess)+": CC*bias-invalid vs corr-cue, "+newline+"r = " + num2str(round(corr_cc_biasInv_cue.r,3))+", p = " + num2str(round(corr_cc_biasInv_cue.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("bias"); ylabel("DBM 'r' ");

set(gcf,'position',[50,50,1800,400])
%linkaxes([ax1 ax2 ax3 ax4],'xy')
saveas(gcf,plot_savepath+type+"_cc_biasInv_rt.png")
