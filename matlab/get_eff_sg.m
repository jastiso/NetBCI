function [ E_high, E_low, E_other] = get_eff_sg( subset, coeff)
% Get distribution of energy for all subjgraph categories,
%   Detailed explanation goes here


% which one is hhighest for behavior
b_exp = subset(:,end);
[~,b_idx] = max(b_exp);
[~,nb_idx] = min(b_exp);
nSG = size(subset,1);

% get matrix
bSG = get_sg_matrix(102, subset(b_idx,1:end-1));
nbSG = get_sg_matrix(102, subset(nb_idx, 1:end-1));

% get efficiency
E_high = efficiency_wei(bSG);
E_low = efficiency_wei(nbSG);
tmp_E = 0;
cnt = 0;
for n = 1:nSG-2
    cnt = cnt + 1;
    curr_SG = get_sg_matrix(102, subset(n,1:end-1));
    E = efficiency_wei(curr_SG);
    if n ~= b_idx && n ~= nb_idx
        tmp_E = tmp_E + E;
    end
end
E_other = tmp_E/(nSG-2);
end
