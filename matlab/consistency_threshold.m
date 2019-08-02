%% Consistency Threshold



%% Global Variables

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/param/';
save_dir = [top_dir, 'GroupAvg/wpli/analysis/'];

subjs = [1:20];
nSubj = numel(subjs);
nNode = 102;
nEdge = (nNode^2-nNode)/2;
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'}];
sensors = [{'grad'}];
% threshold for edge visualization
thr = 1;
consensus = zeros(nNode,nNode, nSubj, numel(bands));
load([save_dir, 'noise_sg.mat']);
noise_idx = noise_sg;
img_dir = [top_dir, 'GroupAvg/wpli/images/'];

% initialize
high_states = zeros(nNode, nSubj, numel(bands));
low_states = zeros(nNode, nSubj, numel(bands));

%% Loop through data

for i = subjs
    subj = sprintf('%03d', i);
    for j = 1:numel(sensors)
        sens = sensors{j};
        load([top_dir, 'idx_', sens, '.mat'])
        for k = 1:numel(bands)
            f = bands{k};
            img_dir = [top_dir, 'NMF/', subj, '/', sens, '/images/'];
            % make directories
            if ~exist(img_dir, 'dir')
                mkdir(img_dir)
            end
            save_dir = [top_dir, 'NMF/', subj, '/', sens, '/analysis/wpli/'];
            % make directories
            if ~exist(save_dir, 'dir')
                mkdir(save_dir)
            end
            
            % get subgraph data
            if strcmp(data_dir(end-5:end-1),'param')
                subset = readNPY([data_dir, subj, '/', sens, 'wpli_', f, '_subset.npy']);
                coeff = readNPY([data_dir, subj, '/',sens, 'wpli_', f, '_coeff.npy']);
            else
                subset = readNPY([data_dir, subj, '/', sens, '/wpli_', f, '_subset.npy']);
                coeff = readNPY([data_dir, subj, '/',sens, '/wpli_', f, '_coeff.npy']);
            end
            
           
            % get subgraph with highest expresssion in behavior
            % remove noise SG
            idx = noise_sg{k,i};
            coeff = coeff(~idx,:);
            subset = subset(~idx,:);
            [~,bsg] = max(subset(:,end));
            
            % sizes
            nSG = size(subset,1);
            nSens = 102;
            
            for n = 1:nSG
                % get expressin into matrix
                curr_SG = subset(n,1:end-1);
                node_exp = get_sg_matrix(nSens, curr_SG);
                if n == bsg
                    consensus(:,:,i,k) = node_exp;
                end
            end
        end
    end
end


%% plot consensus - consistency (but without distance)

if ~exist(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/images/'], 'dir')
    mkdir(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/images/'])
end

thrs = linspace(0.01, 1, 10);
consistent_edges = zeros(numel(thrs), numel(bands));
for k = 1:numel(thrs)
    thrsh_mats = consensus;
    % threshold
    for i = subjs
        for j = 1:numel(bands)
            thrsh_mats(:,:,i,j) = thresh_mat(consensus(:,:,i,j),thrs(k));
        end
    end
    G = squeeze(sum(thrsh_mats > 0, 3));
    for m = 1:numel(bands)
        consistent_edges(k,m) = mean(sum(G(:,:,m)))/2;
    end
end
% plot
figure(1); clf
plot(thrs,consistent_edges, 'linewidth', 2)
legend(bands)
saveas(gca, [img_dir, 'consensus_thresh.png'], 'png')
