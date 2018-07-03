%% Get expression stats
% find edges with most similar expression across graphs via
% correlations



addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

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

%% Loop through data

for i = 2:nSubj
    subj = sprintf('%03d', i);
    for j = 1:numel(sensors)
        sens = sensors{j};
        load([top_dir, 'idx_', sens, '.mat'])
        for k = 1:numel(bands)
            f = bands{k};

            % get subgraph data
            coeff = readNPY([data_dir, subj, '/', f, '_coeff.npy']);
            subset = readNPY([data_dir, subj, '/', f, '_subset.npy']);
                
            b_exp = coeff(:,end);
            edge_similarity(:,i-1,k,j) = corr(b_exp, coeff(:,1:end-1));
        end
    end
end


% plot
thr = 500;
sim_edges = zeros(thr, nSubj-1, numel(bands));
for i = 1:numel(bands)
    [plot_sim, idx] = sort(edge_similarity(:,:,i), 'descend');
    sim_edges(:,:,i) = idx(1:thr,:);
    plot(plot_sim, 'linewidth', 2)
    xlabel('Edge')
    xlim([0,size(plot_sim,1)])
    ylabel('Pearsons Corr (r)')
    saveas(gca, [img_dir, bands{i}, 'corr_edges.png'], 'png')
end

% get similar edges
sim_reduced = cell(1,numel(bands));
for i = 1:numel(bands)
   unqA = unique(sim_edges(:,:,i));
    sim_reduced{i} = unqA(all(any(bsxfun(@eq,sim_edges(:,:,i),permute(unqA,[2,3,1])),1),2));
end
