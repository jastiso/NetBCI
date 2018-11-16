%% Find cross subject structure
% this is an exploratory script, trying to see if it is worth pursueing
% cross subject cluster of subgraphs, rather than grouping them by rank


addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir = '/Users/stiso/Documents/R/NetBCI/data/wpli/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/param/';
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
bands = [{'alpha'}, {'beta'}, {'low_gamma'}];
sensors = [{'grad'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;

% load behavior
load([top_dir, 'Behavior/stats'])


%% Make correlation matrix for every band

load([save_dir, 'noise_sg.mat']);
noise_idx = noise_sg;

corr_mat = cell(3,1);
loading = cell(numel(bands),1);
for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir_s = [R_dir, sens, '/'];
    
    for k = 1:numel(bands)
        f = bands{k};
        
        % collect all subgraphs
        all_sg = [];
        for i = Subj
            s_idx = find(i == Subj);
            subj = sprintf('%03d', i);
            
            % get subgraph data
            if strcmp(data_dir(end-5:end-1),'param')
                subset = readNPY([data_dir, subj, '/', sens, 'wpli_', f, '_subset.npy']);
                coeff = readNPY([data_dir, subj, '/',sens, 'wpli_', f, '_coeff.npy']);
            else
                subset = readNPY([data_dir, subj, '/', sens, '/wpli_', f, '_subset.npy']);
                coeff = readNPY([data_dir, subj, '/',sens, '/wpli_', f, '_coeff.npy']);
            end
            
            % remove noise SG
            idx = noise_sg{k,i};
            coeff = coeff(~idx,:);
            subset = subset(~idx,:);
            
            all_sg = [all_sg; subset(:,1:(end-1))];
            loading{k} = [loading{k}; subset(:,end)];
        end
        corr_mat{k} = corr(all_sg');
        corr_mat{k} = corr_mat{k} - eye(size(corr_mat{k},1));
        figure(1); clf
        imagesc(corr_mat{k}); colorbar
        saveas(gca, [img_dir, f, '_corr_sub.png'], 'png')
    end
end


%% Get clustering

ci = cell(1,numel(bands));
Q = zeros(1,numel(bands));
for k = 1:numel(bands)
    f = bands{k};
    [a,b] = modularity_und(corr_mat{k},1);
    ci{k} = a;
    Q(k) = b;
end

%plot according to cluster
for k = 1:numel(bands)
    curr = corr_mat{k};
    [~,idx] = sort(ci{k});
    figure(k); clf
    imagesc(curr(idx,idx)); colorbar;
    title(['Q = ', num2str(Q(k))])
    saveas(gca, [img_dir, bands{k}, '_corr.png'], 'png')
end


%% Compare behavior loading

for k = 1:numel(bands)
    anova1(loading{k}, ci{k})
    saveas(gca, [img_dir, bands{k}, '_bp.png'], 'png')
end

