function [RT_bin,  pred_avg] = bin_pred_vec(pred_vec,RT,bin_type)

%% Group trials into 3 bins based on the sequence of trials in t, t-1, t-2, & t-3

% 1) bins: split into three groups (low, medium, high) either by sorting
% and dividing into tertiles (pred_event, pred_entropy) or based on the
% absolute value of the predicted probability (abs_pred_event)
% 2) bin_type: 
%               a) pred_event = predicted probability of event (p)
%               b) pred_entropy = p*(1-p)
%               c) abs_pred_event = predicted probability of event; bins
%                  based on absolute thresholds (0---1/3---2/3---1)


%% Code

N = length(pred_vec);

bin_low = 1:round((N/3));
bin_med = round(N/3)+1:round(2*N/3);
bin_high = round(2*N/3)+1:N;

if(bin_type == "pred_event")
    
    [sort_pred,idx_sort_pred] = sort(pred_vec);
    
    RT_temp = RT(idx_sort_pred);
    RT_bin{1} = RT_temp(bin_low);
    RT_bin{2} = RT_temp(bin_med);
    RT_bin{3} = RT_temp(bin_high);
    
    pred_avg(1) = mean(sort_pred(bin_low));
    pred_avg(2) = mean(sort_pred(bin_med));
    pred_avg(3) = mean(sort_pred(bin_high));
    
elseif(bin_type == "pred_entropy")
    
    pred_entropy = pred_vec.*(1-pred_vec);
    [sort_pred_entropy,idx_sort_pred_entropy] = sort(pred_entropy);
    
    RT_temp = RT(idx_sort_pred_entropy);
    RT_bin{1} = RT_temp(bin_low);
    RT_bin{2} = RT_temp(bin_med);
    RT_bin{3} = RT_temp(bin_high);
    
    pred_avg(1) = mean(sort_pred_entropy(bin_low));
    pred_avg(2) = mean(sort_pred_entropy(bin_med));
    pred_avg(3) = mean(sort_pred_entropy(bin_high));
    
elseif(bin_type == "abs_pred_event")
    
    RT_bin{1} = RT(pred_vec<=(1/3));
    RT_bin{2} = RT(pred_vec>(1/3) & pred_vec<=(2/3));
    RT_bin{3} = RT(pred_vec>(2/3));
   
    pred_avg(1) = mean(pred_vec(pred_vec<=(1/3)));
    pred_avg(2) = mean(pred_vec(pred_vec>(1/3) & pred_vec<=(2/3)));
    pred_avg(3) = mean(pred_vec(pred_vec>(2/3)));
end



end