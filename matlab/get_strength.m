function [ St_high, St_low, St_other] = get_strength( subset, coeff)
% Get distribution of energy for all subjgraph categories,
%   Detailed explanation goes here


% which one is hhighest for behavior
b_exp = subset(:,end);
[~,bSG] = max(b_exp);
[~,nbSG] = min(b_exp);
nSG = size(subset,1);

% get sparsity
St_high = median(subset(bSG,1:end-1));
St_low = median(subset(nbSG,1:end-1));
tmp_St = 0;
cnt = 0;
for n = 1:nSG-2
    cnt = cnt + 1;
    St = median(subset(n,1:end-1));
    if n ~= bSG && n ~= nbSG
        tmp_St = tmp_St + St;
    end
end
St_other = tmp_St/(nSG-2);
end




