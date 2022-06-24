%% Bin trials based on the sequence of recent trials and plot the mean predicted probability (DBM) obtained from gridbest search(Fig 1 Yu, Cohen paper)
% Looking for evidence of sequential effects in the pred data
% Trials pooled across blocks / subjects


%% code
clear;
addpath('functions/')
addpath('results/')

side = ["Left","Right"];
sess = ["Sham","Stim","Post"];
var_type = "choiceRA";

data_to_load='results/compareGridbest_normDivideByMeanThenLog_'+var_type+'.mat';
load(data_to_load);

beh_metric = "RT"; % "RT", "acc" or "DBM_pred"

lim.RT = [0.5,0.9];
lim.acc = [0.2,1];
lim.DBM_pred = [0,1];


%% 
for iSide = 1:2
    for iBin = 1:16
       num_choice{iBin}='0';
    end    

    for iSess = 1:3
        
        for iSubject=1:length(Trials_Info)
            choice_rep=bin_vec_choice{iSide,iSess}(iSubject,:);
            if beh_metric=='DBM_pred'
                ptraj_choice=pred_vec_choice_gridbest{iSide,iSess}(iSubject,:);
                ptraj_pred_choice = zeros(size(ptraj_choice));
                ptraj_pred_choice(choice_rep==1) = ptraj_choice(choice_rep==1);
                ptraj_pred_choice(choice_rep==0) = 1-ptraj_choice(choice_rep==0);
                beh_choice = 1-ptraj_pred_choice; %for RT
                %beh_valid = ptraj_pred_valid; %for accuracy
                
            elseif beh_metric=='RT'
                rt=rtvals{iSide,iSess}(iSubject,:);
                beh_choice= rt; 
            end
            
            %to_plot=acc_vals{iSide,iSess}(iSubject,:);
            to_plot=rtvals{iSide,iSess}(iSubject,:);
            step=49;
           
            
            for iBlock=1:5
                start_idx=(iBlock-1)*step+1;
                end_idx=iBlock*step;

%                 [beh_bin_valid{iSubject,iBlock}] = bin_sequence(valid_rep(start_idx:end_idx),beh_valid(start_idx:end_idx),rt(start_idx:end_idx));
%                 [beh_bin_invalid{iSubject,iBlock}] = bin_sequence(invalid_rep(start_idx:end_idx),beh_invalid(start_idx:end_idx),rt(start_idx:end_idx));
%                 [beh_bin_nochange{iSubject,iBlock}] = bin_sequence(nochange_rep(start_idx:end_idx),beh_nochange(start_idx:end_idx),rt(start_idx:end_idx));
%                 [beh_bin_cue{iSubject,iBlock}] = bin_sequence(cue_rep(start_idx:end_idx),beh_cue(start_idx:end_idx),rt(start_idx:end_idx));
%                 
        
                [beh_bin_choice{iSubject,iBlock}] = bin_sequence(choice_rep(start_idx:end_idx),beh_choice(start_idx:end_idx),to_plot(start_idx:end_idx));
               
            end
            
        end
        
        beh_choice = cell(1,16);
         
        for iBin = 1:16
            for iSubject = 1:length(Trials_Info)
                for iBlock = 1:5
                    
                    beh_choice{iBin} = [beh_choice{iBin}; beh_bin_choice{iSubject,iBlock}{iBin}'];
                   
                end
            end
            
           
            meanBeh_choice(iBin) = mean(beh_choice{iBin});
            num_choice{iBin} = num2str(str2num(num_choice{iBin})+length(beh_choice{iBin}));
            
        end
        
        
        %% Plotting
        
        x_label_name = {'1111','0111','1011','0011','1101','0101','1001','0001','1110','0110','1010','0010','1100','0100','1000','0000'};
        if beh_metric=='DBM_pred'
            ylim_val=lim.DBM_pred;
        elseif beh_metric=='RT'
            ylim_val=lim.RT;
        end
        
        ylim_val=[-0.3,0.4];
        
        figure(iSide*100)
        hold on;
        if(iSess == 1)
            plot(meanBeh_choice,'LineWidth',1,'Marker','o','Color',[0,0,1]);
        elseif(iSess==2)
            plot(meanBeh_choice,'LineWidth',1,'Marker','o','Color',[0,0,0]);
        else
            plot(meanBeh_choice,'LineWidth',1,'Marker','o','Color',[1,0,0]);
        end
        ylim(ylim_val)
        ylabel("RT")
        xticklabels(x_label_name);
        xtickangle(65); xticks([1:16]);
        legend({'Sham','Stim','Post'});
        
        set(gcf,'position',[50,100,2000,750])
        
        if iSess==1
            ax1=gca;
            ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
            set(ax2, 'XAxisLocation', 'top');
            set(ax2, 'XLim', get(ax1, 'XLim'));
            set(ax2, 'XTick', get(ax1, 'XTick'));
            set(ax2, 'XTickLabel', num_choice);
            xlabel(beh_metric+" "+side(iSide)+": choice RA");
            xtickangle(65)
        end
        
               
    end
    saveas(gcf,"results/gridbest_seq_effects_normDivideByMeanThenLog_"+beh_metric+"_"+var_type+"_"+side(iSide)+".png");
end

function adjust_plot(lim,label,num_trials,x_label_name)
ylabel("RT")
xticklabels(x_label_name);
xtickangle(65); xticks([1:16]);
legend({'Sham','Stim','Post'});
ax1=gca;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top');
set(ax2, 'XLim', get(ax1, 'XLim'));
set(ax2, 'XTick', get(ax1, 'XTick'));
set(ax2, 'XTickLabel', num_trials);
xlabel(label);
set(ax2, 'YTickLabel', []);
set(ax2, 'YTick', []);
xtickangle(65)
set(gcf,'position',[50,100,2000,750])
ylim(lim)

end


