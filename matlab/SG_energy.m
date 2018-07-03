%% Test the energy of different subgraphs


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
E_high = zeros(nSubj-1,numel(bands));
E_low = zeros(nSubj-1,numel(bands));
E_perf = zeros(nSubj-1, numel(bands));
E_other = zeros(nSubj - 1, numel(bands));
E_exp_corr = [];
high_corr = [];
band_order = [];
subj_order = [];
cnt = 0;
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
            
            
            % get energy
            E_high(i-1,k) = sum(subset(bSG,:).^2);
            E_low(i-1,k) = sum(subset(nbSG,:).^2);
            high_corr = [high_corr; E_low(i-1,k), coeff(nbSG, end)];
            E_perf(i-1,k) = sum(curr_perf.^2)';
            tmp_E = 0;
            for n = 1:nSG-2
                cnt = cnt + 1;
                E = sum(subset(n,:).^2);
                E_exp_corr = [E_exp_corr; E, coeff(n,end)];
                band_order(cnt,1) = k;
                subj_order(cnt,1) = i;
                if n ~= bSG && n ~= nbSG
                    tmp_E = tmp_E + E;
                end
            end
            E_other(i-1,k) = tmp_E/(nSG-2);
        end
    end
end

% plot
for i = 1:numel(bands)
    figure(1); clf
   boxplot([E_high(:,i), E_low(:,i), E_perf(:,i), E_other(:,i)], [{'High'}, {'Low'}, {'Behavior'}, {'Other'}])
   saveas(gca, [save_dir, bands{i}, '_energy.png'], 'png')
end

corrplot(E_exp_corr, 'Var', {'E', 'BE'}, 'type', 'Spearman', 'testR', 'on')
saveas(gca, [save_dir, 'E_BE_corr_all.png'], 'png')
corrplot(high_corr, 'Var', {'E', 'BE'}, 'type', 'Spearman', 'testR', 'on')
saveas(gca, [save_dir, 'E_BE_corr_low.png'], 'png')


% save data
save([R_dir, 'E.mat'], 'E_high', 'E_low', 'E_other');
save([R_dir, 'E_corr.mat'], 'E_exp_corr', 'band_order', 'subj_order');


%% Repeat for Baseline

% Loop through data
E_high = zeros(nSubj-1,numel(bands));
E_low = zeros(nSubj-1,numel(bands));
E_perf = zeros(nSubj-1, numel(bands));
E_other = zeros(nSubj - 1, numel(bands));
E_exp_corr = [];
high_corr = [];
band_order = [];
subj_order = [];
cnt = 0;
for i = 2:nSubj
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
            
            
            % get energy
            E_high(i-1,k) = sum(subset(bSG,:).^2);
            E_low(i-1,k) = sum(subset(nbSG,:).^2);
            high_corr = [high_corr; E_low(i-1,k), coeff(nbSG, end)];
            E_perf(i-1,k) = sum(curr_perf.^2)';
            tmp_E = 0;
            for n = 1:nSG-2
                cnt = cnt + 1;
                E = sum(subset(n,:).^2);
                E_exp_corr = [E_exp_corr; E, coeff(n,end)];
                band_order(cnt,1) = k;
                subj_order(cnt,1) = i;
                if n ~= bSG && n ~= nbSG
                    tmp_E = tmp_E + E;
                end
            end
            E_other(i-1,k) = tmp_E/(nSG-2);
        end
    end
end

% plot
for i = 1:numel(bands)
    figure(1); clf
   boxplot([E_high(:,i), E_low(:,i), E_perf(:,i), E_other(:,i)], [{'High'}, {'Low'}, {'Behavior'}, {'Other'}])
   saveas(gca, [save_dir, bands{i}, '_energy_bl.png'], 'png')
end

%corrplot(E_exp_corr, 'Var', {'E', 'BE'}, 'type', 'Spearman', 'testR', 'on')
%saveas(gca, [save_dir, 'E_BE_corr_all_bl.png'], 'png')
%corrplot(high_corr, 'Var', {'E', 'BE'}, 'type', 'Spearman', 'testR', 'on')
%saveas(gca, [save_dir, 'E_BE_corr_low_bl.png'], 'png')


% save data
save([R_dir, 'E_bl.mat'], 'E_high', 'E_low', 'E_other');
save([R_dir, 'E_corr_bl.mat'], 'E_exp_corr', 'band_order', 'subj_order');