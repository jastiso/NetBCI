function [ H_high, H_low, H_perf, H_other, H_exp_corr, band_order, subj_order ] = get_entropy( subset, coeff, perf, s, b)
% Get distribution of energy for all subjgraph categories,
%   Detailed explanation goes here

H_exp_corr = [];

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
H_high = wentropy(coeff(bSG,:), 'shannon');
H_low = wentropy(coeff(nbSG,:), 'shannon');
high_corr = [H_low, subset(nbSG, end)];
H_perf = sum(curr_perf.^2)';
tmp_H = 0;
cnt = 0;
for n = 1:nSG-2
    cnt = cnt + 1;
    H = wentropy(coeff(n,:), 'shannon');
    H_exp_corr = [H_exp_corr; H, subset(n,end)];
    band_order(cnt,1) = b;
    subj_order(cnt,1) = s;
    if n ~= bSG && n ~= nbSG
        tmp_H = tmp_H + H;
    end
end
H_other = tmp_H/(nSG-2);
end




