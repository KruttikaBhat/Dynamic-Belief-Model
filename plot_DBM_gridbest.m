%% Plots the correlation between RT and DBM predictions across subjects.
% load data to be plotted

%% code

% Sort subjects based on the correlation coefficient (avg across
% left/right) for each variable/event:
clear;
addpath('functions/')
addpath('results/')

var_type = "respRA1";
data_to_load='results/Accuracy_compareGridbest_'+var_type+'.mat';
load(data_to_load);

[~,sort_valid] = sort(rvals_valid{1,1}(:,4)+rvals_valid{2,1}(:,4));
[~,sort_invalid] = sort(rvals_invalid{1,1}(:,4)+rvals_invalid{2,1}(:,4));
[~,sort_nochange] = sort(rvals_nochange{1,1}(:,4)+rvals_nochange{2,1}(:,4));
[~,sort_cue] = sort(rvals_cue{1,1}(:,4)+rvals_cue{2,1}(:,4));

%ylim_val=[-0.8,0.1];
ylim_val=[-1,1];

for iSide = 1:2
    for iSess = 1:3
        figure(iSide*100);  hold on;
        
        subplot(2,2,1); hold on;
        if(iSess == 1)
            plot(rvals_valid{iSide,iSess}(sort_valid,4),'LineWidth',1,'Marker','o','Color',[0,0,1]);
        elseif(iSess==2)
            plot(rvals_valid{iSide,iSess}(sort_valid,4),'LineWidth',1,'Marker','o','Color',[0,0,0]);
        else
            plot(rvals_valid{iSide,iSess}(sort_valid,4),'LineWidth',1,'Marker','o','Color',[1,0,0]);
        end
        title(side(iSide)+": ValidTrial");
        xlabel('Subject'); ylabel('Correlation(RT, pred)');
        ylim(ylim_val)
        
        subplot(2,2,2); hold on;
        if(iSess == 1)
            plot(rvals_invalid{iSide,iSess}(sort_invalid,4),'LineWidth',1,'Marker','o','Color',[0,0,1]);
        elseif(iSess==2)
            plot(rvals_invalid{iSide,iSess}(sort_invalid,4),'LineWidth',1,'Marker','o','Color',[0,0,0]);
        else
            plot(rvals_invalid{iSide,iSess}(sort_invalid,4),'LineWidth',1,'Marker','o','Color',[1,0,0]);
        end
        title(side(iSide)+": InvalidTrial");
        xlabel('Subject'); ylabel('Correlation(RT, pred)');
        ylim(ylim_val)
        
        subplot(2,2,3); hold on;
        if(iSess == 1)
            plot(rvals_nochange{iSide,iSess}(sort_nochange,4),'LineWidth',1,'Marker','o','Color',[0,0,1]);
        elseif(iSess==2)
            plot(rvals_nochange{iSide,iSess}(sort_nochange,4),'LineWidth',1,'Marker','o','Color',[0,0,0]);
        else
            plot(rvals_nochange{iSide,iSess}(sort_nochange,4),'LineWidth',1,'Marker','o','Color',[1,0,0]);
        end
        title(side(iSide)+": NoChangeTrial");
        xlabel('Subject'); ylabel('Correlation(RT, pred)');
        ylim(ylim_val)
        
        subplot(2,2,4); hold on;
        if(iSess == 1)
            plot(rvals_cue{iSide,iSess}(sort_cue,4),'LineWidth',1,'Marker','o','Color',[0,0,1]);
        elseif(iSess==2)
            plot(rvals_cue{iSide,iSess}(sort_cue,4),'LineWidth',1,'Marker','o','Color',[0,0,0]);
        else
            plot(rvals_cue{iSide,iSess}(sort_cue,4),'LineWidth',1,'Marker','o','Color',[1,0,0]);
        end
        title(side(iSide)+": Cue");
        xlabel('Subject'); ylabel('Correlation(RT, pred)');
        ylim(ylim_val)
         
    end
    legend({'Sham','Stim','Post'},'Location','southeast');
    set(gcf,'position',[50,50,1000,800])
    
    saveas(gcf,"results/Accuracy_corr_compareGridbest_"+var_type+"_"+side(iSide)+".png")
    
end
