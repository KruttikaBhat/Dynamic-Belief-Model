clear;clc;

% fname='Ankita_data.json';
% numsubjects=26;
% numsess=6;
% numblocks=5;

%fname='Sanjna_data.json';
%numsubjects=28;
%numsess=2;
%numblocks=5;

%fname='Varsha_data_oddeven.json';
%numsubjects=22;
%numsess=2;
%numblocks=5;

fname='Varsha_data_10Blocks.json';
numsubjects=22;
numsess=1;
numblocks=10;

fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
json_val = jsondecode(str);

%data(no.trials,no.columns,no.blocks)
%sno,validity,cueside,angle,lc,rc,filler_col,filler_col,choice,rt

for subject=1:numsubjects
    for sess_idx=1:numsess
        for block=1:numblocks
            idx=find(json_val.subject_num==subject-1 & json_val.session_idx==sess_idx-1 ...
                & json_val.block_num==block-1);
            all_data(:,1,block,sess_idx,subject)=json_val.trial_num(idx);
            all_data(:,2,block,sess_idx,subject)=json_val.cue_validity(idx);
            all_data(:,3,block,sess_idx,subject)=json_val.cue_location(idx);
            all_data(:,4,block,sess_idx,subject)=json_val.stimuli(idx);
            all_data(:,5,block,sess_idx,subject)=json_val.left_change(idx);
            all_data(:,6,block,sess_idx,subject)=json_val.right_change(idx);
            all_data(:,7,block,sess_idx,subject)=zeros(length(idx),1);
            all_data(:,8,block,sess_idx,subject)=zeros(length(idx),1);
            all_data(:,9,block,sess_idx,subject)=json_val.choice(idx);
            all_data(:,10,block,sess_idx,subject)=json_val.rt(idx);
            
        end
        all_angles=unique(all_data(:,4,:,sess_idx,subject));
        angles(:,sess_idx,subject)=all_angles(all_angles~=0);
        numangles=length(angles(:,sess_idx,subject));
        [contable(:,:,:,sess_idx,subject), contable_lc, contable_rc] = CreateConTable_hemifield(all_data(:,:,:,sess_idx,subject), angles(:,sess_idx,subject));
        
        %zero correction
        temp_con=contable(:,:,:,sess_idx,subject);
        temp_con(temp_con==0)=0.5/numangles;
        contable(:,:,:,sess_idx,subject)=temp_con;
        
        [theta_est, theta_err, LLF, ctable_mod, ctable_fit] = mADC_model_fit(contable(:,:,:,sess_idx,subject));
        d_val(subject,sess_idx,:)=theta_est(1:2:numangles*2);
        d_inv(subject,sess_idx,:)=theta_est(2:2:numangles*2);
        c_val(subject,sess_idx,:)=theta_est(2*numangles+1);
        c_inv(subject,sess_idx,:)=theta_est(2*numangles+2);
        
        cc_val_init(subject,sess_idx,:) = c_val(subject,sess_idx,:) - d_val(subject,sess_idx,:)/2;
        cc_inv_init(subject,sess_idx,:) = c_inv(subject,sess_idx,:) - d_inv(subject,sess_idx,:)/2;
        
        cc_val(subject,sess_idx,:)=mean(cc_val_init(subject,sess_idx,:),3);
        cc_inv(subject,sess_idx,:)=mean(cc_inv_init(subject,sess_idx,:),3);
        
        d_val_avg(subject,sess_idx,:)=mean(d_val(subject,sess_idx,:),3);
        d_inv_avg(subject,sess_idx,:)=mean(d_inv(subject,sess_idx,:),3);

        fprintf("%d %d",subject,sess_idx);
  
    end
    
end

save('Varsha_10Blocks_perf_metrics');
