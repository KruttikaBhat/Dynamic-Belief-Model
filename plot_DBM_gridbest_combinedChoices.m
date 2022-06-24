%% Plots the correlation between RT and DBM predictions across subjects.
% load data to be plotted

%% code

% Sort subjects based on the correlation coefficient (avg across
% left/right) for each variable/event:
clear;
addpath('functions/')
addpath('results/')

var_type = "choiceRA";
data_to_load='results/compareGridbest_normDivideByMeanThenLog_'+var_type+'.mat';
load(data_to_load);

[~,sort_choice] = sort(rvals_choice{1,1}(:,4)+rvals_choice{2,1}(:,4));

%ylim_val=[-0.8,0.1];
ylim_val=[-1,1];

for iSide = 1:2
    for iSess = 1:3
        figure(iSide*100);  hold on;
        
        if(iSess == 1)
            plot(rvals_choice{iSide,iSess}(sort_choice,4),'LineWidth',1,'Marker','o','Color',[0,0,1]);
        elseif(iSess==2)
            plot(rvals_choice{iSide,iSess}(sort_choice,4),'LineWidth',1,'Marker','o','Color',[0,0,0]);
        else
            plot(rvals_choice{iSide,iSess}(sort_choice,4),'LineWidth',1,'Marker','o','Color',[1,0,0]);
        end
        title(side(iSide)+": ChoiceRepetitionTrial");
        xlabel('Subject'); ylabel('Correlation(RT, pred)');
        ylim(ylim_val)
        
         
    end
    legend({'Sham','Stim','Post'},'Location','southeast');
    set(gcf,'position',[50,50,1000,800])
    
    saveas(gcf,"results/RT_corr_normDivideByMeanThenLog_compareGridbest_"+var_type+"_"+side(iSide)+".png")
    
end
