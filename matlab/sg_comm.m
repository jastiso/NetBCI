%% Get communication values for different subgraphs
% we can opporationalize this easily as efficiency and communicability. I
% also want to try navigation metrics, search info, and possibly
% broadcasting?


% @author JStiso jeni.stiso@gmail.com

% Change Log
% June 19 2018, created script

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/BCT/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir = '/Users/stiso/Documents/R/NetBCI/data/gc/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/gc/';
save_dir = [top_dir, 'GroupAvg/gc/analysis/'];
img_dir = [top_dir, 'GroupAvg/gc/images/'];
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
    
    G_high = zeros(nSubj,numel(bands));
    G_low = zeros(nSubj,numel(bands));
    G_other = zeros(nSubj, numel(bands));
    G_exp_corr = [];
    E_high = zeros(nSubj,numel(bands));
    E_low = zeros(nSubj,numel(bands));
    E_other = zeros(nSubj, numel(bands));
    E_exp_corr = [];
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
            subset = readNPY([data_dir, subj, '/', sens, '/gc_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/gc_', f, '_coeff.npy']);
            b_exp = subset(:,end);
            [~,bSG] = max(b_exp);
            [~,nbSG] = min(b_exp);
            nSG = size(subset,1);
            
            % communicability
            [ high, low, other] = get_G_sg(subset, coeff);
            G_high(s_idx,k) = high; G_low(s_idx,k) = low; G_other(s_idx,k) = other;
            clear high low e_perf other exp_corr
           
            % efficiency
            [ high, low, other] = get_eff_sg(subset, coeff);
            E_high(s_idx,k) = high; E_low(s_idx,k) = low; E_other(s_idx,k) = other;
            clear high low e_perf other exp_corr

        end
    end
    
    
    % plot
    % G
    for i = 1:numel(bands)
        figure(1); clf
        boxplot([G_high(:,i), G_low(:,i), G_other(:,i)], [{'High'}, {'low'}, {'Other'}])
        saveas(gca, [save_dir, bands{i}, '_G.png'], 'png')
    end
    
    % save data
    save([R_dir_s, 'G.mat'], 'G_high', 'G_low', 'G_other');
   
    
    % efficiency
    for n = 1:numel(bands)
        figure(1); clf
        boxplot([E_high(:,n), E_low(:,n), E_other(:,n)], [{'High'}, {'Low'}, {'Other'}])
        saveas(gca, [save_dir, bands{n}, '_eff.png'], 'png')
    end

    % save data
    save([R_dir_s, 'Eff.mat'], 'E_high', 'E_low', 'E_other');
    
end
