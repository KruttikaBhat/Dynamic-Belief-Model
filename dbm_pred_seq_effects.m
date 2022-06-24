%% Bin trials based on the sequence of recent trials and plot the mean predicted probability (DBM) obtained from gridbest search(Fig 1 Yu, Cohen paper)
% Looking for evidence of sequential effects in the pred data
% Trials pooled across blocks / subjects


%% code
clear;
addpath('functions/')
addpath('results/')

side = ["Left","Right"];
sess = ["Sham","Stim","Post"];
var_type = "respRA2";

data_to_load='results/compareGridbest_Accuracy_'+var_type+'.mat';
load(data_to_load);

beh_metric = "DBM_pred"; % "RT" or "DBM_pred"

lim.RT = [0.5,0.9];
lim.DBM_pred = [0,1];


%% 
for iSide = 1:2
    
    for iBin = 1:16
       num_valid{iBin}='0';
       num_invalid{iBin}='0';
       num_nochange{iBin}='0';
       num_cue{iBin}='0';
    end    

    for iSess = 1:3
        
        for iSubject=1:length(Trials_Info)
            valid_rep=bin_vec_valid{iSide,iSess}(iSubject,:);
            invalid_rep=bin_vec_invalid{iSide,iSess}(iSubject,:);
            nochange_rep=bin_vec_nochange{iSide,iSess}(iSubject,:);
            cue_rep=bin_vec_cue{iSide,iSess}(iSubject,:);
            
            if beh_metric=='DBM_pred'
                ptraj_valid=pred_vec_valid_gridbest{iSide,iSess}(iSubject,:);
                ptraj_pred_valid = zeros(size(ptraj_valid));
                ptraj_pred_valid(valid_rep==1) = ptraj_valid(valid_rep==1);
                ptraj_pred_valid(valid_rep==0) = 1-ptraj_valid(valid_rep==0);
                beh_valid = 1-ptraj_pred_valid; %for RT
                %beh_valid = ptraj_pred_valid; %for accuracy
                    
                ptraj_invalid=pred_vec_invalid_gridbest{iSide,iSess}(iSubject,:);
                ptraj_pred_invalid = zeros(size(ptraj_invalid));
                ptraj_pred_invalid(invalid_rep==1) = ptraj_invalid(invalid_rep==1);
                ptraj_pred_invalid(invalid_rep==0) = 1-ptraj_invalid(invalid_rep==0);
                beh_invalid = 1-ptraj_pred_invalid;
                %beh_invalid = ptraj_pred_invalid;
                
                ptraj_nochange = pred_vec_nochange_gridbest{iSide,iSess}(iSubject,:);
                ptraj_pred_nochange = zeros(size(ptraj_nochange));
                ptraj_pred_nochange(nochange_rep==1) = ptraj_nochange(nochange_rep==1);
                ptraj_pred_nochange(nochange_rep==0) = 1-ptraj_nochange(nochange_rep==0);
                beh_nochange = 1-ptraj_pred_nochange;
                %beh_nochange = ptraj_pred_nochange;
                    
                ptraj_cue = pred_vec_cue_gridbest{iSide,iSess}(iSubject,:);
                ptraj_pred_cue = zeros(size(ptraj_cue));
                ptraj_pred_cue(cue_rep==1) = ptraj_cue(cue_rep==1);
                ptraj_pred_cue(cue_rep==0) = 1-ptraj_cue(cue_rep==0);
                beh_cue = 1-ptraj_pred_cue;
                %beh_cue = ptraj_pred_cue;
                
            elseif beh_metric=='RT'
                rt=rtvals{iSide,iSess}(iSubject,:);
                beh_valid= rt; beh_invalid= rt; beh_nochange = rt; beh_cue = rt;
            end
            
            %to_plot=acc_vals{iSide,iSess}(iSubject,:);
            to_plot=rtvals{iSide,iSess}(iSubject,:);
           
            if(var_type=="respRA1" || var_type=="respRA2")
                step=49;
            else
                step=50;
            end
            
            for iBlock=1:5
                start_idx=(iBlock-1)*step+1;
                end_idx=iBlock*step;

%                 [beh_bin_valid{iSubject,iBlock}] = bin_sequence(valid_rep(start_idx:end_idx),beh_valid(start_idx:end_idx),rt(start_idx:end_idx));
%                 [beh_bin_invalid{iSubject,iBlock}] = bin_sequence(invalid_rep(start_idx:end_idx),beh_invalid(start_idx:end_idx),rt(start_idx:end_idx));
%                 [beh_bin_nochange{iSubject,iBlock}] = bin_sequence(nochange_rep(start_idx:end_idx),beh_nochange(start_idx:end_idx),rt(start_idx:end_idx));
%                 [beh_bin_cue{iSubject,iBlock}] = bin_sequence(cue_rep(start_idx:end_idx),beh_cue(start_idx:end_idx),rt(start_idx:end_idx));
%                 
        
                [beh_bin_valid{iSubject,iBlock}] = bin_sequence(valid_rep(start_idx:end_idx),beh_valid(start_idx:end_idx),to_plot(start_idx:end_idx));
                [beh_bin_invalid{iSubject,iBlock}] = bin_sequence(invalid_rep(start_idx:end_idx),beh_invalid(start_idx:end_idx),to_plot(start_idx:end_idx));
                [beh_bin_nochange{iSubject,iBlock}] = bin_sequence(nochange_rep(start_idx:end_idx),beh_nochange(start_idx:end_idx),to_plot(start_idx:end_idx));
                [beh_bin_cue{iSubject,iBlock}] = bin_sequence(cue_rep(start_idx:end_idx),beh_cue(start_idx:end_idx),to_plot(start_idx:end_idx));

            end
            
        end
        
        beh_valid = cell(1,16);
        beh_invalid = cell(1,16);
        beh_nochange = cell(1,16);
        beh_cue = cell(1,16);
         
        for iBin = 1:16
            for iSubject = 1:length(Trials_Info)
                for iBlock = 1:5
                    
                    beh_valid{iBin} = [beh_valid{iBin}; beh_bin_valid{iSubject,iBlock}{iBin}'];
                    beh_invalid{iBin} = [beh_invalid{iBin}; beh_bin_invalid{iSubject,iBlock}{iBin}'];
                    beh_nochange{iBin} = [beh_nochange{iBin}; beh_bin_nochange{iSubject,iBlock}{iBin}'];
                    beh_cue{iBin} = [beh_cue{iBin}; beh_bin_cue{iSubject,iBlock}{iBin}'];
                    
                end
            end
            
           
            meanBeh_valid(iBin) = mean(beh_valid{iBin});
            meanBeh_invalid(iBin) = mean(beh_invalid{iBin});
            meanBeh_nochange(iBin) = mean(beh_nochange{iBin});
            meanBeh_cue(iBin) = mean(beh_cue{iBin});
           
            num_valid{iBin} = num2str(str2num(num_valid{iBin})+length(beh_valid{iBin}));
            num_invalid{iBin} = num2str(str2num(num_invalid{iBin})+length(beh_invalid{iBin}));
            num_nochange{iBin} = num2str(str2num(num_nochange{iBin})+length(beh_nochange{iBin}));
            num_cue{iBin} = num2str(str2num(num_cue{iBin})+length(beh_cue{iBin}));
          
        end
        
        
        %% Plotting
        
        x_label_name = {'1111','0111','1011','0011','1101','0101','1001','0001','1110','0110','1010','0010','1100','0100','1000','0000'};
        if beh_metric=='DBM_pred'
            ylim_val=lim.DBM_pred;
        elseif beh_metric=='RT'
            ylim_val=lim.RT;
        end
        
        ylim_val=[0,1];
        
        figure(iSide*100)
        subplot(2,2,1)
        hold on;
        if(iSess == 1)
            plot(meanBeh_valid,'LineWidth',1,'Marker','o','Color',[0,0,1]);
        elseif(iSess==2)
            plot(meanBeh_valid,'LineWidth',1,'Marker','o','Color',[0,0,0]);
        else
            plot(meanBeh_valid,'LineWidth',1,'Marker','o','Color',[1,0,0]);
        end
        ylim(ylim_val)
        adjust_plot(lim.(beh_metric),beh_metric+" "+side(iSide)+": valid",num_valid,x_label_name);
        
        subplot(2,2,2)
        hold on;
        if(iSess == 1)
            plot(meanBeh_invalid,'LineWidth',1,'Marker','o','Color',[0,0,1]);
        elseif(iSess==2)
            plot(meanBeh_invalid,'LineWidth',1,'Marker','o','Color',[0,0,0]);
        else
            plot(meanBeh_invalid,'LineWidth',1,'Marker','o','Color',[1,0,0]);
        end
        ylim(ylim_val)
        adjust_plot(lim.(beh_metric),beh_metric+" "+side(iSide)+": invalid",num_invalid,x_label_name);

        subplot(2,2,3)
        hold on;
        if(iSess == 1)
            plot(meanBeh_nochange,'LineWidth',1,'Marker','o','Color',[0,0,1]);
        elseif(iSess==2)
            plot(meanBeh_nochange,'LineWidth',1,'Marker','o','Color',[0,0,0]);
        else
            plot(meanBeh_nochange,'LineWidth',1,'Marker','o','Color',[1,0,0]);
        end
        ylim(ylim_val)
        adjust_plot(lim.(beh_metric),beh_metric+" "+side(iSide)+": nochange",num_nochange,x_label_name);

        subplot(2,2,4)
        hold on;
        if(iSess == 1)
            plot(meanBeh_cue,'LineWidth',1,'Marker','o','Color',[0,0,1]);
        elseif(iSess==2)
            plot(meanBeh_cue,'LineWidth',1,'Marker','o','Color',[0,0,0]);
        else
            plot(meanBeh_cue,'LineWidth',1,'Marker','o','Color',[1,0,0]);
        end
        ylim(ylim_val)
        adjust_plot(lim.(beh_metric),beh_metric+" "+side(iSide)+": cue",num_cue,x_label_name);

               
    end
    saveas(gcf,"results/gridbest_seq_effects_Accuracy_"+beh_metric+"_"+var_type+"_"+side(iSide)+".png");
end

function adjust_plot(lim,label,num_trials,x_label_name)
ylabel("accuracy")
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


