function [ G_high, G_low, G_other] = get_G_sg( subset, coeff)
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

% get weighted communicability
Db = diag(sum(bSG));
Dnb = diag(sum(nbSG));
G_high = mean(mean(expm(Db^(-1/2)*bSG*Db^(-1/2))));
G_low = mean(mean(expm(Dnb^(-1/2)*nbSG*Dnb^(-1/2))));
tmp_G = 0;
cnt = 0;
for n = 1:nSG-2
    cnt = cnt + 1;
    curr_SG = get_sg_matrix(102, subset(n,1:end-1));
    D = diag(sum(curr_SG));
    G = sum(mean(mean(D^(-1/2)*curr_SG*D^(-1/2))));
    if n ~= b_idx && n ~= nb_idx
        tmp_G = tmp_G + G;
    end
end
G_other = tmp_G/(nSG-2);
end

