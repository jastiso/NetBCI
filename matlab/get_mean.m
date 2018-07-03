function [ M_high, M_low, M_perf, M_other, M_exp_corr, band_order, subj_order ] = get_mean( subset, coeff, perf, s, b)
% Get distribution of energy for all subjgraph categories,
%   Detailed explanation goes here

M_exp_corr = [];

% which one is hhighest for behavior
b_exp = subset(:,end);
[~,bSG] = max(b_exp);
[~,nbSG] = min(b_exp);
nSG = size(subset,1);

% scale perf to what was used in NMF
curr_perf = perf(:,s);
curr_perf = curr_perf/mean(curr_perf);
curr_perf = curr_perf.*mean(mean(coeff));

% get energy
M_high = mean(coeff(bSG,:));
M_low = mean(coeff(nbSG,:));
high_corr = [M_low, subset(nbSG, end)];
M_perf = sum(curr_perf.^2)';
tmp_S = 0;
cnt = 0;
for n = 1:nSG-2
    cnt = cnt + 1;
    S = mean(coeff(n,:));
    M_exp_corr = [M_exp_corr; S, subset(n,end)];
    band_order(cnt,1) = b;
    subj_order(cnt,1) = s;
    if n ~= bSG && n ~= nbSG
        tmp_S = tmp_S + S;
    end
end
M_other = tmp_S/(nSG-2);
end




