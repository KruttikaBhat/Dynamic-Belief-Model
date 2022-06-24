function [RT_bin, cnt] = bin_sequence(rep_vec,RT,RT_act)

%% Group trials into 16 bins based on the sequence of trials in t, t-1, t-2, & t-3

% See Fig 1A (Yu, Cohen 2008) - same bins/analysis used here
% Format: Bin_Num = (t-3)(t-2)(t-1)(t)

% Bins:

bin_seq{1} = [1,1,1,1];
bin_seq{2} = [0,1,1,1];
bin_seq{3} = [1,0,1,1];
bin_seq{4} = [0,0,1,1];

bin_seq{5} = [1,1,0,1];
bin_seq{6} = [0,1,0,1];
bin_seq{7} = [1,0,0,1];
bin_seq{8} = [0,0,0,1];

bin_seq{9}  = [1,1,1,0];
bin_seq{10} = [0,1,1,0];
bin_seq{11} = [1,0,1,0];
bin_seq{12} = [0,0,1,0];

bin_seq{13} = [1,1,0,0];
bin_seq{14} = [0,1,0,0];
bin_seq{15} = [1,0,0,0];
bin_seq{16} = [0,0,0,0];

cnt = ones(1,16);
RT_bin = cell(1,16);

for t = 4:length(rep_vec)
    
    trial_seq = rep_vec(t-3:t);
    
    for iBin = 1:16
        if(trial_seq == bin_seq{iBin} & RT_act(t)~=0)
        %if(trial_seq == bin_seq{iBin})
            if(RT_act(t)==0)
                disp(RT(t))
            end    
            RT_bin{iBin}(cnt(iBin)) = RT(t);
            cnt(iBin) = cnt(iBin)+1;
        end
    end
    
end

end