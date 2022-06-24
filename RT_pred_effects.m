%% Bin trials into low, medium, high groups based on predicted probabilities (from DBM) and check if RT is different across these bins

% Similar to raw_seq_effects code. Check that for more detailed comments.

% 1) bins: split into three groups (low, medium, high) either by sorting
% and dividing into tertiles (pred_event, pred_entropy) or based on the
% absolute value of the predicted probability (abs_pred_event)
% 2) bin_type: 
%               a) pred_event = predicted probability of event (p)
%               b) pred_entropy = p*(1-p)
%               c) abs_pred_event = predicted probability of event; bins
%                  based on absolute thresholds (0---1/3---2/3---1)

%% Code

load ('data/tACS_40Hz_woETrej.mat')

%load('results/compareGridbest_xyz.mat') % Load the DBM results file containing the predicted probability vectors.

side = ["Left","Right"];
sess = ["Sham","Stim","Post"];
var_type = "resp";
bin_type = "pred_event"; 


for iSide = 1:2
    
    if(iSide == 1)
        Trials_Info = leftPPC.Trials_Info;
    else
        Trials_Info = rightPPC.Trials_Info;
    end
    
    for iSess = 1:3
        
        
        for iSubject=1:length(Trials_Info)
            data_sub = Trials_Info{iSubject};
            
            if(iSess==1)
                all_blocks = data_sub.Sham;
            elseif(iSess==2)
                all_blocks = data_sub.Stim;
            elseif(iSess==3)
                all_blocks = data_sub.Post;
            end
            
            for iBlock = 1:5
                ablock = all_blocks(:,:,iBlock);
                if(var_type == "actual")
                    isvalid_cue = ablock(:, 2) == 2;
                    isleft_cue = ablock(:,3) == -1;
                    ischange_trial = ablock(:,5) == 1 | ablock(:,6) == 1;
                    
                elseif(var_type == "resp")
                    isvalid_cue = ablock(:, 3) == ablock(:,9) & ablock(:,9)~=0; % cue side == resp side; resp side ~= no change
                    isleft_cue = ablock(:,3) == -1;
                    ischange_trial = ablock(:,9) ~= 0; %% resp ~= no change
                    
                end
                
                for iTrial = 2:size(ablock,1)
                    %                     isvalid_cue_rep(iTrial) = isvalid_cue(iTrial) == isvalid_cue(iTrial-1);
                    %                     isleft_cue_rep(iTrial) = isleft_cue(iTrial) == isleft_cue(iTrial-1);
                    %                     ischange_trial_rep(iTrial) = ischange_trial(iTrial) == ischange_trial(iTrial-1);
                    
                    isvalid_cue_rep(iTrial) = isvalid_cue(iTrial) ;
                    isleft_cue_rep(iTrial) = isleft_cue(iTrial) ;
                    ischange_trial_rep(iTrial) = ischange_trial(iTrial) ;
                end
                rt{iSubject}((iBlock-1)*50+1:(iBlock-1)*50+50) = ablock(:,10)./max(ablock(:,10));
                
                
            end
            [RT_bin_isvalid{iSubject}, avg_pred_isvalid(iSubject,:)] = bin_pred_vec(pred_vec_isvalid{iSide,iSess}(iSubject,:),rt{iSubject},bin_type);
            [RT_bin_isleft{iSubject}, avg_pred_isleft(iSubject,:)] = bin_pred_vec(pred_vec_isleft_gridbest{iSide,iSess}(iSubject,:),rt{iSubject},bin_type);
            [RT_bin_ischange{iSubject}, avg_pred_ischange(iSubject,:)] = bin_pred_vec(pred_vec_ischange_gridbest{iSide,iSess}(iSubject,:),rt{iSubject},bin_type);
           
            
        end
        
        RT_isvalid = cell(1,3);
        RT_isleft = cell(1,3);
        RT_ischange = cell(1,3);
        
        for iBin = 1:3
            for iSubject = 1:length(Trials_Info)
                
                
                RT_isvalid{iBin} = [RT_isvalid{iBin}; RT_bin_isvalid{iSubject}{iBin}'];
                RT_isleft{iBin} = [RT_isleft{iBin}; RT_bin_isleft{iSubject}{iBin}'];
                RT_ischange{iBin} = [RT_ischange{iBin}; RT_bin_ischange{iSubject}{iBin}'];
                
                
                
            end
            MRT_isvalid(iBin) = median(RT_isvalid{iBin});
            MRT_isleft(iBin) = median(RT_isleft{iBin});
            MRT_ischange(iBin) = median(RT_ischange{iBin});
            num_isvalid{iBin} = num2str(length(RT_isvalid{iBin}));
            num_isleft{iBin} = num2str(length(RT_isleft{iBin}));
            num_ischange{iBin} = num2str(length(RT_ischange{iBin}));
            
        end
        
        
        x_label_name = {'Low','Medium','High'};
        
        figure(iSide*100)
        subplot(1,3,1)
        hold on;
        if(iSess == 1)
            plot(MRT_isvalid,'LineWidth',3,'Color',[0,0,1])
            
        elseif(iSess==2)
            plot(MRT_isvalid,'LineStyle','--','LineWidth',3,'Color',[0.3010, 0.7450, 0.9330])
            
        else
            plot(MRT_isvalid,'Color','r')
            
        end
        ylabel("RT")
        xticklabels(x_label_name);
        xtickangle(65); xticks([1:10]);
        ylim([0.44,0.66])
        
        ax1=gca;
        ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
        set(ax2, 'XAxisLocation', 'top');
        set(ax2, 'XLim', get(ax1, 'XLim'));
        set(ax2, 'XTick', get(ax1, 'XTick'));
        set(ax2, 'XTickLabel', num_isvalid)
        set(ax2, 'YTickLabel', []);
        set(ax2, 'YTick', []);
        
        xlabel(side(iSide)+": isvalid");
        xtickangle(65)
        
        subplot(1,3,2)
        hold on;
        
        if(iSess == 1)
            plot(MRT_isleft,'LineWidth',3,'Color',[0,0,1])
            
        elseif(iSess==2)
            plot(MRT_isleft,'LineStyle','--','LineWidth',3,'Color',[0.3010, 0.7450, 0.9330])
            
        else
            plot(MRT_isleft,'Color','r')
            
        end
        ylabel("RT")
        xticklabels(x_label_name);
        xtickangle(65); xticks([1:3]);
        ylim([0.44,0.66])
        
        ax1=gca;
        ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
        set(ax2, 'XAxisLocation', 'top');
        set(ax2, 'XLim', get(ax1, 'XLim'));
        set(ax2, 'XTick', get(ax1, 'XTick'));
        set(ax2, 'XTickLabel', num_isleft)
        xlabel(side(iSide)+": isleft");
        set(ax2, 'YTickLabel', []);
        set(ax2, 'YTick', []);
        
        xtickangle(65)
        
        subplot(1,3,3)
        hold on;
        
        if(iSess == 1)
            plot(MRT_ischange,'LineWidth',3,'Color',[0,0,1])
            
        elseif(iSess==2)
            plot(MRT_ischange,'LineStyle','--','LineWidth',3,'Color',[0.3010, 0.7450, 0.9330])
            
        else
            plot(MRT_ischange,'Color','r')
            
        end
        ylabel("RT")
        xticklabels(x_label_name);
        xtickangle(65); xticks([1:3]);
        ylim([0.44,0.66])
        legend({'Sham','Stim','Post'});
        
        ax1=gca;
        ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
        set(ax2, 'XAxisLocation', 'top');
        set(ax2, 'XLim', get(ax1, 'XLim'));
        set(ax2, 'XTick', get(ax1, 'XTick'));
        set(ax2, 'XTickLabel', num_ischange)
        xlabel(side(iSide)+": ischange");
        set(ax2, 'YTickLabel', []);
        set(ax2, 'YTick', []);
        
        xtickangle(65)
        set(gcf,'position',[50,100,1600,800])
        
    end
    
    saveas(gcf,"results/RT_pred_effects_"+var_type+"_"+side(iSide)+".png");
end
