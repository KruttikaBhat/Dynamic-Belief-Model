function rt_means=get_subject_means(leftPPC,rightPPC)


for iSubject=1:26
    rt_vals=[];
    for iSide = 1:2
    
        if(iSide == 1)
            Trials_Info = leftPPC.Trials_Info;
        else
            Trials_Info = rightPPC.Trials_Info;
        end
    
        for iSess = 1:3
        
            data_sub = Trials_Info{iSubject};
            
            if(iSess==1)
                all_blocks = data_sub.Sham;
            elseif(iSess==2)
                all_blocks = data_sub.Stim;
            elseif(iSess==3)
                all_blocks = data_sub.Post;
            end
            for block=1:size(all_blocks,3)
                for trl=1:size(all_blocks,1)
                    if all_blocks(trl,10)~=0
                        rt_vals=[rt_vals;all_blocks(trl,10)];
                    end
                end
            end
        end
    end
    rt_means{iSubject}=mean(rt_vals);
end
end


