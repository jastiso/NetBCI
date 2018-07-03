%% Test the skew of different coeffs


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
S_high = zeros(nSubj-1,numel(bands));
S_low = zeros(nSubj-1,numel(bands));
S_other = zeros(nSubj - 1, numel(bands));
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
            S_high(i-1,k) = skewness(coeff(bSG,1:end-1));
            S_low(i-1,k) = skewness(coeff(nbSG,1:end-1));
            tmp_S = 0;
            for n = 1:nSG-2
                cnt = cnt + 1;
                S = skewness(coeff(n,1:end-1));
                band_order(cnt,1) = k;
                subj_order(cnt,1) = i;
                if n ~= bSG && n ~= nbSG
                    tmp_S = tmp_S + S;
                end
            end
            S_other(i-1,k) = tmp_S/(nSG-2);
        end
    end
end

% plot
for i = 1:numel(bands)
    figure(1); clf
   boxplot([S_high(:,i), S_low(:,i), S_other(:,i)], [{'High'}, {'Low'}, {'Other'}])
   saveas(gca, [save_dir, bands{i}, '_skew_coeff.png'], 'png')
end


% save data
save([R_dir, 'S_coeff.mat'], 'S_high', 'S_low', 'S_other');

%% Now for baseline

% Loop through data
S_high = zeros(nSubj-1,numel(bands));
S_low = zeros(nSubj-1,numel(bands));
S_other = zeros(nSubj - 1, numel(bands));
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
            S_high(i-1,k) = skewness(coeff(bSG,1:end-1));
            S_low(i-1,k) = skewness(coeff(nbSG,1:end-1));
            tmp_S = 0;
            for n = 1:nSG-2
                cnt = cnt + 1;
                S = skewness(coeff(n,1:end-1));
                band_order(cnt,1) = k;
                subj_order(cnt,1) = i;
                if n ~= bSG && n ~= nbSG
                    tmp_S = tmp_S + S;
                end
            end
            S_other(i-1,k) = tmp_S/(nSG-2);
        end
    end
end

% plot
for i = 1:numel(bands)
    figure(1); clf
   boxplot([S_high(:,i), S_low(:,i), S_other(:,i)], [{'High'}, {'Low'}, {'Other'}])
   saveas(gca, [save_dir, bands{i}, '_skew_coeff_bl.png'], 'png')
end


% save data
save([R_dir, 'S_coeff_bl.mat'], 'S_high', 'S_low', 'S_other');

