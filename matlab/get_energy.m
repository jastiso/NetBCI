function [ E_high, E_low, E_perf, E_other, E_exp_corr, band_order, subj_order ] = get_energy( subset, coeff, perf, s, b)
% Get distribution of energy for all subjgraph categories,

% Change Log
% July 10 - updated this to run statistics for all but the highest...since
% the highest tends to look exaclty like the behavior

E_exp_corr = [];

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
E_high = sum(coeff(bSG,:).^2);
E_low = sum(coeff(nbSG,:).^2);
high_corr = [E_low, subset(nbSG, end)];
E_perf = sum(curr_perf.^2)';
tmp_E = 0;
cnt = 0;
for n = 1:nSG-2
    cnt = cnt + 1;
    E = sum(coeff(n,:).^2);
    band_order(cnt,1) = b;
    subj_order(cnt,1) = s;
    if n ~= bSG
        E_exp_corr = [E_exp_corr; E, subset(n,end)];
        if n ~= nbSG
            tmp_E = tmp_E + E;
        end
    end
end
E_other = tmp_E/(nSG-2);
end

