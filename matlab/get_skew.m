function [ S_high, S_low, S_perf, S_other, S_exp_corr, band_order, subj_order ] = get_skew( subset, coeff, perf, s, b)
% Get distribution of energy for all subjgraph categories,

% Change Log
% July 10 - updated this to run statistics for all but the highest...since
% the highest tends to look exaclty like the behavior

S_exp_corr = [];

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
S_high = skewness(coeff(bSG,:));
S_low = skewness(coeff(nbSG,:));
high_corr = [S_low, subset(nbSG, end)];
S_perf = sum(curr_perf.^2)';
tmp_S = 0;
cnt = 0;
for n = 1:nSG-2
    cnt = cnt + 1;
    S = skewness(coeff(n,:));
    band_order(cnt,1) = b;
    subj_order(cnt,1) = s;
    S_exp_corr = [S_exp_corr; S, subset(n,end)];
    if n ~= nbSG && n ~= bSG
        tmp_S = tmp_S + S;
    end
    
end
S_other = tmp_S/(nSG-2);
end




