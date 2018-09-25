%% run sub graph statistics for wpli



addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir = '/Users/stiso/Documents/R/NetBCI/data/wpli/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/';
save_dir = [top_dir, 'GroupAvg/wpli/analysis/'];
img_dir = [top_dir, 'GroupAvg/wpli/images/'];
% make directories
if ~exist(img_dir, 'dir')
    mkdir(img_dir)
end
% make directories
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
% make directories
if ~exist(R_dir, 'dir')
    mkdir(R_dir)
end

Subj = [1:20];
nSubj = numel(Subj);
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'}, {'beta'}, {'low_gamma'}, {'gamma'}];
sensors = [{'grad'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;
edge_similarity = zeros(nEdges, nSubj-1, numel(bands), numel(sensors));
% load behavior
load([top_dir, 'Behavior/perf'])

%% Loop through data
for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir_s = [R_dir, sens, '/'];
    
    St_high = zeros(nSubj,numel(bands));
    St_low = zeros(nSubj,numel(bands));
    St_other = zeros(nSubj, numel(bands));
    St_exp_corr = [];
    S_high = zeros(nSubj-1,numel(bands));
    S_low = zeros(nSubj-1,numel(bands));
    S_other = zeros(nSubj - 1, numel(bands));
    S_exp_corr = [];
    band_order = [];
    subj_order = [];
    cnt = 0;
    for i = Subj
        s_idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        load([top_dir, 'idx_', sens, '.mat'])
        for k = 1:numel(bands)
            f = bands{k};
            
            % get subgraph data
            subset = readNPY([data_dir, subj, '/', sens, '/wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/wpli_', f, '_coeff.npy']);
            b_exp = coeff(:,end);
            [~,bSG] = max(b_exp);
            [~,nbSG] = min(b_exp);
            nSG = size(subset,1);
            
            % energy
            [ high, low, other] = get_strength(subset, coeff);
            St_high(s_idx,k) = high; St_low(s_idx,k) = low; St_other(s_idx,k) = other;
            clear high low e_perf other exp_corr
           
            %skew
            [ high, low, other] = get_skew_sg(subset, coeff);
            S_high(s_idx,k) = high; S_low(s_idx,k) = low; S_other(s_idx,k) = other;
            clear high low e_perf other exp_corr

        end
    end
    
    
    % plot
    % strength
    for i = 1:numel(bands)
        figure(1); clf
        boxplot([St_high(:,i), St_low(:,i), St_other(:,i)], [{'High'}, {'Behavior'}, {'Other'}])
        saveas(gca, [save_dir, bands{i}, '_sparsity.png'], 'png')
    end
    
    % save data
    save([R_dir_s, 'St.mat'], 'St_high', 'St_low', 'St_other');
   
    
    % skew
    for n = 1:numel(bands)
        figure(1); clf
        boxplot([S_high(:,n), S_low(:,n), S_other(:,n)], [{'High'}, {'Low'}, {'Other'}])
        saveas(gca, [save_dir, bands{n}, '_skew_sg.png'], 'png')
    end

    % save data
    save([R_dir_s, 'Ssg.mat'], 'S_high', 'S_low', 'S_other');
    
end
