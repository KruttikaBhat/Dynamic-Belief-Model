% Same as corr_pred_sens_bias, but taking delta-sens (sens(cued) - sens(uncued)) or delta-bias


%% code
clear;
type="respRA2";
sess=["Sham"];
plot_savepath="results/corr_pred_normDivideByMeanThenLog_";

%get Ankita's data
load("../DBM_fit_Ankita_data/results/compareGridbest_normDivideByMeanThenLog_"+type+".mat") % Load the DBM results file to be plotted
load("../DBM_fit_Ankita_data/data/Ankita_perf_metrics.mat")
iSess=1;
dataset=1;

for iSide = 1:2
    for iSubject = 1:26            
        sensDelta{dataset}((iSide-1)*26+iSubject) = mean(d_val(iSubject,(iSide-1)*3+iSess,:))-mean(d_inv(iSubject,(iSide-1)*3+iSess,:));
        cc_biasDelta{dataset}((iSide-1)*26+iSubject) = cc_val(iSubject,(iSide-1)*3+iSess)-cc_inv(iSubject,(iSide-1)*3+iSess);

        predCorr_val{dataset}((iSide-1)*26+iSubject) = rvals_valid{iSide,iSess}(iSubject,4);
        predCorr_inval{dataset}((iSide-1)*26+iSubject) = rvals_invalid{iSide,iSess}(iSubject,4);
        predCorr_nochange{dataset}((iSide-1)*26+iSubject) = rvals_nochange{iSide,iSess}(iSubject,4);
        predCorr_cue{dataset}((iSide-1)*26+iSubject) = rvals_cue{iSide,iSess}(iSubject,4);

    end
end
  
%get Sanjna's data
load("../DBM_fit_Sanjna_data/results/compareGridbest_normDivideByMeanThenLog_"+type+".mat") % Load the DBM results file to be plotted
load("../DBM_fit_Sanjna_data/data/Sanjna_perf_metrics.mat")

iSess=1;
dataset=2;

for iSubject = 1:28
    sensDelta{dataset}(iSubject) = mean(d_val(iSubject,iSess,:))-mean(d_inv(iSubject,iSess,:));
    cc_biasDelta{dataset}(iSubject) = cc_val(iSubject,iSess)-cc_inv(iSubject,iSess);
    
    predCorr_val{dataset}(iSubject) = rvals_valid{iSess}(iSubject,4);
    predCorr_inval{dataset}(iSubject) = rvals_invalid{iSess}(iSubject,4);
    predCorr_nochange{dataset}(iSubject) = rvals_nochange{iSess}(iSubject,4);
    predCorr_cue{dataset}(iSubject) = rvals_cue{iSess}(iSubject,4);
end


%get Varsha's data
load("../DBM_fit_Varsha_data/results/10Blocks/compareGridbest_normDivideByMeanThenLog_"+type+".mat") % Load the DBM results file to be plotted
load("../DBM_fit_Varsha_data/data/Varsha_10Blocks_perf_metrics.mat")

iSess=1;
dataset=3;

for iSubject = 1:22
    sensDelta{dataset}(iSubject) = mean(d_val(iSubject,iSess,:))-mean(d_inv(iSubject,iSess,:));
    cc_biasDelta{dataset}(iSubject) = cc_val(iSubject,iSess)-cc_inv(iSubject,iSess);

    predCorr_val{dataset}(iSubject) = rvals_valid{iSess}(iSubject,4);
    predCorr_inval{dataset}(iSubject) = rvals_invalid{iSess}(iSubject,4);
    predCorr_nochange{dataset}(iSubject) = rvals_nochange{iSess}(iSubject,4);
    predCorr_cue{dataset}(iSubject) = rvals_cue{iSess}(iSubject,4);
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
title(sess(iSess)+": delta-sens vs corr-valid, "+newline+"r = " + num2str(round(corr_sensDelta_valid.r,3))+", p = " + num2str(round(corr_sensDelta_valid.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("delta sens"); ylabel("DBM 'r' ");

subplot(1,4,2); hold on;
scatter(sensDelta{1}(:),predCorr_inval{1}(:),'k','o')
scatter(sensDelta{2}(:),predCorr_inval{2}(:),'k','s')
scatter(sensDelta{3}(:),predCorr_inval{3}(:),'k','^')
title(sess(iSess)+": delta-sens vs corr-invalid, "+newline+"r = " + num2str(round(corr_sensDelta_invalid.r,3))+", p = " + num2str(round(corr_sensDelta_invalid.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("delta sens"); ylabel("DBM 'r' ");

subplot(1,4,3); hold on;
scatter(sensDelta{1}(:),predCorr_nochange{1}(:),'r','o')
scatter(sensDelta{2}(:),predCorr_nochange{2}(:),'r','s')
scatter(sensDelta{3}(:),predCorr_nochange{3}(:),'r','^')
title(sess(iSess)+": delta-sens vs corr-nochange, "+newline+"r = " + num2str(round(corr_sensDelta_nochange.r,3))+", p = " + num2str(round(corr_sensDelta_nochange.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("delta sens"); ylabel("DBM 'r' ");    

subplot(1,4,4); hold on;
scatter(sensDelta{1}(:),predCorr_cue{1}(:),'g','o')
scatter(sensDelta{2}(:),predCorr_cue{2}(:),'g','s')
scatter(sensDelta{3}(:),predCorr_cue{3}(:),'g','^')
title(sess(iSess)+": delta-sens vs corr-cue, "+newline+"r = " + num2str(round(corr_sensDelta_cue.r,3))+", p = " + num2str(round(corr_sensDelta_cue.p,3)));
xlim(xlim_val_d); 
ylim(ylim_val);
xlabel("delta sens"); ylabel("DBM 'r' ");


set(gcf,'position',[50,50,1800,400])
%linkaxes([ax1 ax2 ax3 ax4],'xy')
saveas(gcf,plot_savepath+type+"_delta-sens_rt.png")    

figure(2)
hold on;

subplot(1,4,1); hold on;
scatter(cc_biasDelta{1}(:),predCorr_val{1}(:),'b','o')
scatter(cc_biasDelta{2}(:),predCorr_val{2}(:),'b','s')
scatter(cc_biasDelta{3}(:),predCorr_val{3}(:),'b','^')
title(sess(iSess)+": delta-CC bias vs corr-valid, "+newline+"r = " + num2str(round(corr_cc_biasDelta_valid.r,3))+", p = " + num2str(round(corr_cc_biasDelta_valid.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("delta bias"); ylabel("DBM 'r' ");

subplot(1,4,2); hold on;
scatter(cc_biasDelta{1}(:),predCorr_inval{1}(:),'k','o')
scatter(cc_biasDelta{2}(:),predCorr_inval{2}(:),'k','s')
scatter(cc_biasDelta{3}(:),predCorr_inval{3}(:),'k','^')
title(sess(iSess)+": delta-CC bias vs corr-invalid, "+newline+"r = " + num2str(round(corr_cc_biasDelta_invalid.r,3))+", p = " + num2str(round(corr_cc_biasDelta_invalid.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("delta bias"); ylabel("DBM 'r' ");

subplot(1,4,3); hold on;
scatter(cc_biasDelta{1}(:),predCorr_nochange{1}(:),'r','o')
scatter(cc_biasDelta{2}(:),predCorr_nochange{2}(:),'r','s')
scatter(cc_biasDelta{3}(:),predCorr_nochange{3}(:),'r','^')
title(sess(iSess)+": delta-CC bias vs corr-nochange, "+newline+"r = " + num2str(round(corr_cc_biasDelta_nochange.r,3))+", p = " + num2str(round(corr_cc_biasDelta_nochange.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("delta bias"); ylabel("DBM 'r' ");

subplot(1,4,4); hold on;
scatter(cc_biasDelta{1}(:),predCorr_cue{1}(:),'g','o')
scatter(cc_biasDelta{2}(:),predCorr_cue{2}(:),'g','s')
scatter(cc_biasDelta{3}(:),predCorr_cue{3}(:),'g','^')
title(sess(iSess)+": delta-CC bias vs corr-cue, "+newline+"r = " + num2str(round(corr_cc_biasDelta_cue.r,3))+", p = " + num2str(round(corr_cc_biasDelta_cue.p,3)));
xlim(xlim_val_cc); 
ylim(ylim_val);
xlabel("delta bias"); ylabel("DBM 'r' ");

set(gcf,'position',[50,50,1800,400])
%linkaxes([ax1 ax2 ax3 ax4],'xy')
saveas(gcf,plot_savepath+type+"_delta-cc_bias_rt.png")
