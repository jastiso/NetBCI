%% Get the frequency contribution to each subgraph

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

subjs = [2:5];
f = 'all';
nBand = 4; %theta alpha beta gamma
colors = [{'Pink'}, {'LightGreen'}, {'PaleTurquoise'}, {'Plum'}];
sensors = [{'mag'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;

% load behavior
load([top_dir, 'Behavior/perf'])


%% Loop through subj

for i = subjs
    subj = sprintf('%03d', i);
    for j = 1:numel(sensors)
        sens = sensors{j};
        load([top_dir, 'idx_', sens, '.mat'])
        
        % get subgraph data
        subset = readNPY([data_dir, subj, '/', f, '_subset.npy']);
        coeff = readNPY([data_dir, subj, '/', f, '_coeff.npy']);
        nSG = size(subset,1);
        
        % which one is hhighest for behavior
        b_exp = coeff(:,end);
        band_idx = zeros(size(coeff(:,1:end-1)));
        cnt = 0;
        for k = 1:nEdges:size(band_idx,2)
            cnt = cnt + 1;
            band_idx(:,k:(k+nEdges-1)) = cnt;
        end
        band_idx = [band_idx, zeros(nSG)];
        
        % scale perf to what was used in NMF
        curr_perf = perf(:,i-1);
        curr_perf = curr_perf/mean(curr_perf);
        curr_perf = curr_perf.*mean(mean(subset));
        
        % get distributions for each band
        Ns = zeros(nBand, 50);
        Xs = zeros(nBand, 50);
        for n = 1:nSG
            for k = 1:nBand
                curr_idx = band_idx(n,:) == k;
                tmp = coeff(n,curr_idx);
                %tmp = tmp(tmp ~= 0);
                [N, X] = hist(tmp',50);
                Ns(k,:) = N;
                Xs(k,:) = X;
            end
            % plot
            figure(1); clf; hold on
            for k = 1:nBand
                bar(Xs(k,:), Ns(k,:), 0.8, 'facecolor', rgb(colors{k}),...
                    'edgecolor', rgb(colors{k}),'facealpha', .5);
            end
            title(['BE ', num2str(b_exp(n))])
            xlabel('Coeff'); ylabel('Frequency');
            legend([{'theta'}, {'alpha'}, {'beta'}, {'gamama'}]);
            saveas(gca, [data_dir, subj, '/' f, '_', num2str(n), '_dist.png'], 'png')
        end
        
    end
end