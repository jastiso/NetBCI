%% Plot expression distributions

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

R_dir_o = '/Users/stiso/Documents/R/NetBCI/data/gc/';
top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/gc/';
save_dir = [top_dir, 'GroupAvg/gc/analysis/'];

% make directories
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
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
    load([top_dir, 'idx_', sens, '.mat'])
    R_dir = [R_dir_o, sens, '/'];
    
    S_high = zeros(nSubj-1,numel(bands));
    S_low = zeros(nSubj-1,numel(bands));
    S_perf = zeros(nSubj-1, numel(bands));
    S_other = zeros(nSubj - 1, numel(bands));
    S_exp_corr = [];
    low_corr = [];
    cnt = 0;
    band_order = [];
    subj_order = [];
    for i = Subj
        s_idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        for k = 1:numel(bands)
            f = bands{k};
            
            img_dir = [top_dir, 'NMF/', subj, '/images/', sens, '/'];
            % make directories
            if ~exist(img_dir, 'dir')
                mkdir(img_dir)
            end
            % get subgraph data
            subset = readNPY([data_dir, subj, '/', sens, '/gc_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/', sens, '/gc_', f, '_coeff.npy']);
            nSG = size(subset,1);
            color_sf = 1/(nSG);
            
            % for every band, plot dist of each graph, ranked by BE,
            % highest to lowest - highest exp is red, lowest is blue in
            % plots
            [~, idx] = sort(subset(:,end), 'descend');
            for n = 1:nSG
                eval(['[N', num2str(n), ', X', num2str(n), '] = hist((subset(idx(n),:)),50);'])
                eval(['[cN', num2str(n), ', cX', num2str(n), '] = hist((coeff(idx(n),:)),15);'])
            end
            figure(1); clf
            for n = 1:nSG
                eval(['X = X', num2str(n), ';']);
                eval(['N = N', num2str(n), ';']);
                bar(X1, N, 0.8, 'facecolor', [0+color_sf*n, 0, 1-color_sf*n], 'edgecolor', [0+color_sf*n, 0, 1-color_sf*n],'facealpha', .5); hold on
            end
            title([f, ' SG Exp Histogram'])
            xlabel('(SG weight)')
            hold off
            saveas(gca, [img_dir, f, '_SG_dist.png'], 'png')
            
            figure(2); clf
            for n = 1:nSG
                eval(['cX = cX', num2str(n), ';']);
                eval(['cN = cN', num2str(n), ';']);
                bar(cX1, cN, 0.8, 'facecolor', [0+color_sf*n, 0, 1-color_sf*n], 'edgecolor', [0+color_sf*n, 0, 1-color_sf*n],'facealpha', .5); hold on
            end
            title([f, ' SG Temporal Coeff Histogram'])
            xlabel('temp exp')
            hold off
            saveas(gca, [img_dir, f, '_SG_coeff.png'], 'png')
        end
    end
end