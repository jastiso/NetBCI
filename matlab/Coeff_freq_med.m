%% Get median of each band's distribution

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
bands = [{'theta'},{'alpha'},{'beta'},{'gamma'}];
colors = [{'Pink'}, {'LightGreen'}, {'PaleTurquoise'}, {'Plum'}];
sensors = [{'mag'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;

% load behavior
load([top_dir, 'Behavior/perf'])



%% Loop through data

med_high = zeros(numel(subjs),nBand);
med_low = zeros(numel(subjs),nBand);
med_other = zeros(numel(subjs),nBand);
cnt = 0;
band_order = [];
subj_order = [];

for i = subjs
    subj = sprintf('%03d', i);
    for j = 1:numel(sensors)
        sens = sensors{j};
        load([top_dir, 'idx_', sens, '.mat'])
        
        % get subgraph data
        subset = readNPY([data_dir, subj, '/', f, '_subset.npy']);
        coeff = readNPY([data_dir, subj, '/', f, '_coeff.npy']);
        nSG = size(subset,1);
        
        % band index
        cntb = 0;
        band_idx = zeros(size(coeff(:,1:end-1)));
        for k = 1:nEdges:size(band_idx,2)
            cntb = cntb + 1;
            band_idx(:,k:(k+nEdges-1)) = cntb;
        end
        band_idx = [band_idx, zeros(nSG)];
        
        % which one is hhighest for behavior
        b_exp = coeff(:,end);
        [~,bSG] = max(b_exp);
        [~,nbSG] = min(b_exp);
        
        % scale perf to what was used in NMF
        curr_perf = perf(:,i);
        curr_perf = curr_perf/mean(curr_perf);
        curr_perf = curr_perf.*mean(mean(subset));
        
        
        % get median
        for k = 1:nBand
            
            tmp_med = 0;
            for n = 1:nSG-2
                curr_idx = band_idx(n,:) == k;
                med_high(i,k) = median(coeff(bSG,curr_idx));
                med_low(i,k) = median(coeff(nbSG,curr_idx));
                cnt = cnt + 1;
                med = median(coeff(n,curr_idx));
                band_order(cnt,1) = k;
                subj_order(cnt,1) = i;
                if n ~= bSG && n ~= nbSG
                    tmp_med = tmp_med + med;
                end
            end
            med_other(i,k) = tmp_med/(nSG-2);
        end
    end
end

% plot
for i = 1:nBand
    figure(1); clf
    boxplot([med_high(:,i), med_low(:,i), med_other(:,i)], [{'High'}, {'Low'}, {'Other'}])
    saveas(gca, [save_dir, f, '_med_freq_coeff.png'], 'png')
end


% save data
save([R_dir, 'med_freq_coeff.mat'], 'med_high', 'med_low', 'med_other');

