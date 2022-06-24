function [lr] = get_lr(d,c)

for i = 1:length(d)

lr_temp(i) = exp(d(i)*(c-(0.5*d(i))));

end

lr = 1/(mean(lr_temp));
end