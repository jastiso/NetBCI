function [ Sp_high, Sp_low, Sp_other] = get_skew_sg( subset, coeff)
% Get distribution of energy for all subjgraph categories,
%   Detailed explanation goes here


% which one is hhighest for behavior
b_exp = subset(:,end);
[~,bSG] = max(b_exp);
[~,nbSG] = min(b_exp);
nSG = size(subset,1);

% get sparsity
Sp_high = skewness(subset(bSG,1:end-1));
Sp_low = skewness(subset(nbSG,1:end-1));
tmp_Sp = 0;
cnt = 0;
for n = 1:nSG-2
    cnt = cnt + 1;
    Sp = skewness(subset(n,1:end-1));
    if n ~= bSG && n ~= nbSG
        tmp_Sp = tmp_Sp + Sp;
    end
end
Sp_other = tmp_Sp/(nSG-2);
end





