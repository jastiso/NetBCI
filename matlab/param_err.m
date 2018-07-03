%% Get distributions of parameters and errors

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
nSeed = 100;
nEdges = (nNode^2-nNode)/2;
edge_similarity = zeros(nEdges, nSubj-1, numel(bands), numel(sensors));

%% Loop through data

err_all = zeros(nSubj-1, numel(bands),nSeed);
err_hist = zeros(nSubj-1,numel(bands));
alpha_all = zeros(nSubj-1,numel(bands));
beta_all = zeros(nSubj-1,numel(bands));
n_all = zeros(nSubj-1,numel(bands));

for i = 2:nSubj
    subj = sprintf('%03d', i);
    for j = 1:numel(sensors)
        sens = sensors{j};
        load([top_dir, 'idx_', sens, '.mat'])
        for k = 1:numel(bands)
            f = bands{k};
            
            % get subgraph data
            error = readNPY([data_dir, subj, '/', f, '_err.npy']);
            %readNPY([data_dir, subj, '/', f, '_params.npy']);
            
            % add to vector
            err_hist(i-1,k) = error(end);
            err_all(i-1,k,:) = error;
            
        end
    end
end

%% Plot

% plot convergence
for i = 1:numel(bands)
    figure(1); clf
    plot(squeeze(err_all(:,i,:))', 'linewidth', 2)
    title('Erorr Convergence')
    xlabel('Seeds')
    ylabel('Error')
    saveas(gca, [img_dir, bands{i}, '_err_convergence.png'], 'png')
end

% plot histograms of error, one hist for each band
[N1, X1] = hist(err_all(:,1),15);
[N2, X2] = hist(err_all(:,2),15, X1);
[N3, X3] = hist(err_all(:,3),15, X1);
[N4, X4] = hist(err_all(:,4),15, X1);

figure(1); clf;
bar(X1, N1, 0.8, 'facecolor', rgb('LightSkyBlue'), 'edgecolor', rgb('LightSkyBlue'),'facealpha', .5); hold on
bar(X1, N2, 0.8, 'facecolor', rgb('DeepSkyBlue'), 'edgecolor', rgb('DeepSkyBlue'),'facealpha', .5);
bar(X1, N3, 0.8, 'facecolor', rgb('DodgerBlue'), 'edgecolor', rgb('DodgerBlue'),'facealpha', .5);
bar(X1, N4, 0.8, 'facecolor', rgb('DarkBlue'), 'edgecolor', rgb('DarkBlue'),'facealpha', .5);
title('Final Error Values'); xlabel('Error');
ylabel('Frequency'); legend(bands, 'location', 'NW')
    saveas(gca, [img_dir, 'err_hist.png'], 'png')
