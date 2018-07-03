%% Test the entropy of different subgraphs


addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

R_dir = '/Users/stiso/Documents/R/NetBCI/data/';
top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/';
save_dir = [top_dir, 'GroupAvg/analysis/'];
img_dir = [top_dir, 'GroupAvg/images/'];
% make directories
if ~exist(img_dir, 'dir')
    mkdir(img_dir)
end
% make directories
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

nSubj = 7;
freqs = [3,6;7,14;15,30;40,69];
bands = [{'theta'},{'alpha'},{'beta'},{'gamma'}];
sensors = [{'mag'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;
edge_similarity = zeros(nEdges, nSubj-1, numel(bands), numel(sensors));
% load behavior
load([top_dir, 'Behavior/perf'])

%% Loop through data
H_high = zeros(nSubj-1,numel(bands));
H_low = zeros(nSubj-1,numel(bands));
H_perf = zeros(nSubj-1, numel(bands));
H_other = zeros(nSubj - 1, numel(bands));
H_exp_corr = [];
low_corr = [];
cnt = 0;
band_order = [];
subj_order = [];
for i = 2:nSubj
    subj = sprintf('%03d', i);
    for j = 1:numel(sensors)
        sens = sensors{j};
        load([top_dir, 'idx_', sens, '.mat'])
        for k = 1:numel(bands)
            f = bands{k};

            % get subgraph data
            subset = readNPY([data_dir, subj, '/', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/', f, '_coeff.npy']);
            nSG = size(subset,1);
            
            % which one is hhighest for behavior
            b_exp = coeff(:,end);
            [~,bSG] = max(b_exp);
            [~,nbSG] = min(b_exp);
            
            % scale perf to what was used in NMF
            curr_perf = perf(:,i-1);
            curr_perf = curr_perf/mean(curr_perf);
            curr_perf = curr_perf.*mean(mean(subset));
            
            
            % get shannon entropy
            H_high(i-1,k) = wentropy(subset(bSG,:), 'shannon');
            H_low(i-1,k) = wentropy(subset(nbSG,:), 'shannon');
            low_corr = [low_corr; H_low(i-1,k), coeff(nbSG, end)];
            H_perf(i-1,k) = wentropy(curr_perf', 'shannon');
            tmp_H = 0;
            for n = 1:nSG-2
                cnt = cnt + 1;
                H = wentropy(subset(n,:), 'shannon');
                H_exp_corr = [H_exp_corr; H, coeff(n,end)];
                band_order(cnt,1) = k;
                subj_order(cnt,1) = i;
                if n ~= bSG && n ~= nbSG
                    tmp_H = tmp_H + H;
                end
            end
            H_other(i-1,k) = tmp_H/(nSG-2);
        end
    end
end

% plot
for i = 1:numel(bands)
    figure(1); clf
   boxplot([H_high(:,i), H_low(:,i), H_perf(:,i), H_other(:,i)], [{'High'}, {'Low'}, {'Behavior'}, {'Other'}])
   saveas(gca, [save_dir, bands{i}, '_entropy.png'], 'png')
end

corrplot(H_exp_corr, 'Var', {'H', 'BE'}, 'type', 'Spearman', 'testR', 'on')
saveas(gca, [save_dir, 'H_BE_corr_all.png'], 'png')
corrplot(low_corr, 'Var', {'H', 'BE'}, 'type', 'Spearman', 'testR', 'on')
saveas(gca, [save_dir, 'H_BE_corr_low.png'], 'png')

% save data
save([R_dir, 'H.mat'], 'H_high', 'H_low', 'H_other');
save([R_dir, 'H_corr.mat'], 'H_exp_corr', 'band_order', 'subj_order');


%% Repeat for baseline

% Loop through data
H_high = zeros(nSubj-1,numel(bands));
H_low = zeros(nSubj-1,numel(bands));
H_perf = zeros(nSubj-1, numel(bands));
H_other = zeros(nSubj - 1, numel(bands));
H_exp_corr = [];
low_corr = [];
cnt = 0;
band_order = [];
subj_order = [];
for i = 2:5
    subj = sprintf('%03d', i);
    for j = 1:numel(sensors)
        sens = sensors{j};
        load([top_dir, 'idx_', sens, '.mat'])
        for k = 1:numel(bands)
            f = bands{k};

            % get subgraph data
            subset = readNPY([data_dir, subj, '/', f, '_subset_bl.npy']);
            coeff = readNPY([data_dir, subj, '/', f, '_coeff_bl.npy']);
            nSG = size(subset,1);
            
            % which one is hhighest for behavior
            b_exp = coeff(:,end);
            [~,bSG] = max(b_exp);
            [~,nbSG] = min(b_exp);
            
            % scale perf to what was used in NMF
            curr_perf = perf(:,i-1);
            curr_perf = curr_perf/mean(curr_perf);
            curr_perf = curr_perf.*mean(mean(subset));
            
            
            % get shannon entropy
            H_high(i-1,k) = wentropy(subset(bSG,:), 'shannon');
            H_low(i-1,k) = wentropy(subset(nbSG,:), 'shannon');
            low_corr = [low_corr; H_low(i-1,k), coeff(nbSG, end)];
            H_perf(i-1,k) = wentropy(curr_perf', 'shannon');
            tmp_H = 0;
            for n = 1:nSG-2
                cnt = cnt + 1;
                H = wentropy(subset(n,:), 'shannon');
                H_exp_corr = [H_exp_corr; H, coeff(n,end)];
                band_order(cnt,1) = k;
                subj_order(cnt,1) = i;
                if n ~= bSG && n ~= nbSG
                    tmp_H = tmp_H + H;
                end
            end
            H_other(i-1,k) = tmp_H/(nSG-2);
        end
    end
end

% plot
for i = 1:numel(bands)
    figure(1); clf
   boxplot([H_high(:,i), H_low(:,i), H_perf(:,i), H_other(:,i)], [{'High'}, {'Low'}, {'Behavior'}, {'Other'}])
   saveas(gca, [save_dir, bands{i}, '_entropy_bl.png'], 'png')
end

corrplot(H_exp_corr, 'Var', {'H', 'BE'}, 'type', 'Spearman', 'testR', 'on')
saveas(gca, [save_dir, 'H_BE_corr_all_bl.png'], 'png')
corrplot(low_corr, 'Var', {'H', 'BE'}, 'type', 'Spearman', 'testR', 'on')
saveas(gca, [save_dir, 'H_BE_corr_low_bl.png'], 'png')

% save data
save([R_dir, 'H_bl.mat'], 'H_high', 'H_low', 'H_other');
save([R_dir, 'H_corr_bl.mat'], 'H_exp_corr', 'band_order', 'subj_order');

