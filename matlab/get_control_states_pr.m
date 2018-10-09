%% Get largest eigenvalue state
% get the "easiest to get to" state from the max
% Which channels are used for control: gradiometers:
%GRAD2
%'MEG0242';'MEG0232';'MEG0442';'MEG0432';'MEG0712';'MEG0742';'MEG1842';'MEG1822';'MEG1812';'MEG1622';'MEG1612';
%GRAD3
%'MEG0243';'MEG0233';'MEG0443';'MEG0433';'MEG0713';'MEG0743';'MEG1843';'MEG1823';'MEG1813';'MEG1623';'MEG1613';


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

% things for making control set
grad2 = ['MEG0242';'MEG0232';'MEG0442';'MEG0432';'MEG0712';'MEG0742';'MEG1842';'MEG1822';'MEG1812';'MEG1622';'MEG1612'];
grad3 = ['MEG0243';'MEG0233';'MEG0443';'MEG0433';'MEG0713';'MEG0743';'MEG1843';'MEG1823';'MEG1813';'MEG1623';'MEG1613'];

Subj = [1:20];
nSubj = numel(Subj);
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'}, {'beta'}, {'low_gamma'}, {'gamma'}];
sensors = [{'grad'}];
regions = [{'Left_frontal'}, {'Left_occipital'}, {'Left_parietal'}, {'Left_temporal'}, ...
    {'Right_frontal'}, {'Right_occipital'}, {'Right_parietal'}, {'Right_temporal'}, {'Vertex'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;
load([save_dir, 'pr_noise_sg.mat']);
noise_idx = noise_sg;

% initialize
high_states = zeros(nNode, nSubj, numel(bands));
low_states = zeros(nNode, nSubj, numel(bands));

%% Make control set

load(['/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/2_MEG/1_Signals/2_Segmentation/2_MEG/Session1/1_test01/Seg_MEG_Subj001_Ses1_test01.mat'])
labels = Seg_MEG_Subj001_Ses1_test01.MI.label;
load(['/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/Session1/1_test01/RawData_MEG_Subj001_Ses1_test01.mat'])
mag_idx = strcmp(RawData_MEG_Subj001_Ses1_test01.meg_sensors.chantype,'megmag');
grad_label = labels(~mag_idx);

cnt = 1;
cmb_labels = cell(nNode,1);
for i = 1:2:numel(grad_label)
    cmb_labels{cnt} = [grad_label{i}, '+', grad_label{i+1}(end-3:end)];
    cnt = cnt + 1;
end
lm_label = cell(size(grad2,1),i);
for i = 1:size(grad2,1)
    lm_label{i} = [grad2(i,:), '+', grad3(i,end-3:end)];
end

B = zeros(nNode,1);
for i = 1:size(grad2,1)
    idx = find(contains(cmb_labels,lm_label{i}));
    B(idx) = 1;
end


%% Loop through data

for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir_s = [R_dir, sens, '/'];
    
    for i = Subj
        s_idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        for k = 1:numel(bands)
            f = bands{k};
            
            % get subgraph data
            subset = readNPY([data_dir, subj, '/', sens, 'wpli_pr_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, 'wpli_pr_', f, '_coeff.npy']);
            
            % remove noise SG
            idx = noise_sg{k,i};
            coeff = coeff(~idx,:);
            subset = subset(~idx,:);
            
            b_exp = subset(:,end);
            [~,bSG] = max(b_exp);
            [~,nbSG] = min(b_exp);
            nSG = size(subset,1);
            
            % get gramian
            % scale to be stable
            high_mat = get_sg_matrix(nNode, subset(bSG,:));
            high_mat = high_mat./eigs(high_mat,1) - eye(nNode).*1.001;
            low_mat = get_sg_matrix(nNode, subset(nbSG,:));
            low_mat = low_mat./eigs(low_mat,1) - eye(nNode).*1.001;
            
            % get largest eigenvector of grammian
            [h_vect, h_val] = easy_state(high_mat,diag(B),eye(nNode),zeros(nNode));
            [l_vect, l_val] = easy_state(low_mat,diag(B),eye(nNode),zeros(nNode));
            
            high_states(:,i,k) = h_vect;
            low_states(:,i,k) = l_vect;
        end
    end
end

%% visualize

% figure(1); clf
% imagesc(high_states(:,:,2)); colorbar
%
% figure(2); clf
% imagesc(low_states(:,:,2)); colorbar
%
% figure(3); clf
% imagesc(high_states(:,:,2) - low_states(:,:,2)); colorbar
%
% h_corr = corr(high_states(:,:,2));
% l_corr = corr(low_states(:,:,2));
% figure(1); clf
% imagesc(h_corr); colorbar
%
% figure(2); clf
% imagesc(l_corr); colorbar
%
% figure(1); clf
% imagesc(mean(high_states(:,:,2),2)); colorbar
%
% figure(2); clf
% imagesc(mean(low_states(:,:,2),2)); colorbar

cfg.layout = [top_dir, 'layouts/neuromag306cmb.lay'];


for i = 1:numel(bands)
    plot_data.powspctrm = mean(high_states(:,:,i),2);
    plot_data.label = cmb_labels;
    plot_data.dimord = 'chan_freq';
    plot_data.freq = mean(freqs(i,:));
    plot_data.cfg = [];
    figure(1); clf
    ft_topoplotER(cfg,plot_data); colorbar
    saveas(gca, [img_dir, bands{i}, '_avg_high_state_pr.png'], 'png')
    
    plot_data.powspctrm = mean(low_states(:,:,i),2);
    figure(2); clf
    ft_topoplotER(cfg,plot_data); colorbar
    saveas(gca, [img_dir, bands{i}, '_avg_low_state_pr.png'], 'png')
end

%% Categorize by lobe


h_region = [];
l_region = [];
region_ord = {};
band_ord = {};

cnt = 1;
for i = 1:numel(regions)
    for j = 1:numel(bands)
        load([top_dir, 'montages/', regions{i}, '_idx.mat'])
        
        % get the mean edge between a lobes (total edges/number of edges), divided by
        % the total edges
        
        % stupid frmatting thing for R
        for k = 0:19
            region_ord{cnt+k} = regions{i};
            band_ord{cnt+k} = bands{j};
        end
        h_region = [h_region; mean(high_states(idx,:,j))];
        l_region = [l_region, mean(low_states(idx,:,j))];
        
        cnt = cnt + nSubj;
    end
end

save([R_dir, 'grad/control_state_pr.mat'], 'region_ord', 'band_ord', 'h_region', 'l_region')