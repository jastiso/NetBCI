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
bands = [{'alpha'}, {'beta'}, {'low_gamma'}];
sensors = [{'grad'}];
regions = [{'Left_frontal'}, {'Left_occipital'}, {'Left_parietal'}, {'Left_temporal'}, ...
    {'Right_frontal'}, {'Right_occipital'}, {'Right_parietal'}, {'Right_temporal'}, {'Vertex'}, {'Left_motor'}, {'Right_motor'}];
nNode = 102;
nEdges = (nNode^2-nNode)/2;
load([save_dir, 'pr_noise_sg.mat']);

% initialize
high_states = zeros(nNode, nSubj, numel(bands));
high2_states = zeros(nNode, nSubj, numel(bands));
high3_states = zeros(nNode, nSubj, numel(bands));
high_states_ev2 = zeros(nNode, nSubj, numel(bands));
low_states = zeros(nNode, nSubj, numel(bands));
low_states_ev2 = zeros(nNode, nSubj, numel(bands));
zero_states = zeros(nNode, nSubj, numel(bands));
zero_states_ev2 = zeros(nNode, nSubj, numel(bands));

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
            if strcmp(data_dir(end-5:end-1),'param')
                subset = readNPY([data_dir, subj, '/', sens, 'wpli_pr_', f, '_subset.npy']);
                coeff = readNPY([data_dir, subj, '/',sens, 'wpli_pr_', f, '_coeff.npy']);
            else
                subset = readNPY([data_dir, subj, '/', sens, '/wpli_pr_', f, '_subset.npy']);
                coeff = readNPY([data_dir, subj, '/',sens, '/wpli_pr_', f, '_coeff.npy']);
            end
            
            % remove noise SG
            idx = noise_sg{k,i};
            coeff = coeff(~idx,:);
            subset = subset(~idx,:);
            
            b_exp = subset(:,end);
            [tmp,tmp_idx] = sort(b_exp,'descend');
            bSG = tmp_idx(1);
            bSG2 = tmp_idx(2); % second biggest element
            bSG3 = tmp_idx(3); % third biggest element
            [~, nzSG] = min(nonzeros((b_exp))); % smallest nonzero - this is how we will opporationalize "low"
            
            % if you actualy want to look at 0s
            if sum(b_exp == 0) <= 1
                [~,nbSG] = min(b_exp);
            else
                nbSG = find(b_exp == 0);
            end
            nZero = numel(nbSG);
            
            
            %get gramian
            % scale to be stable
            % 1st
            high_mat = get_sg_matrix(nNode, subset(bSG,:));
            high_mat = high_mat./eigs(high_mat,1) - eye(nNode).*1.001;
            %2nd
            high2_mat = get_sg_matrix(nNode, subset(bSG2,:));
            high2_mat = high2_mat./eigs(high2_mat,1) - eye(nNode).*1.001;
            %3rd
            high3_mat = get_sg_matrix(nNode, subset(bSG3,:));
            high3_mat = high3_mat./eigs(high3_mat,1) - eye(nNode).*1.001;
            %lowest
            low_mat = get_sg_matrix(nNode, subset(nzSG,:));
            low_mat = low_mat./eigs(low_mat,1) - eye(nNode).*1.001;
            zero_mat = zeros(nNode, nNode, nZero);
            for j = 1:nZero
                zero_mat(:,:,j) = get_sg_matrix(nNode, subset(nbSG,:));
                zero_mat(:,:,j) = zero_mat(:,:,j)./eigs(zero_mat(:,:,j),1) - eye(nNode).*1.001;
            end
            
            
            % get largest eigenvector of grammian
            [h_vect, ~, ~, vect2] = easy_state(high_mat,diag(B),eye(nNode),zeros(nNode));
            high_states(:,i,k) = h_vect;
            high_states_ev2(:,i,k) = vect2;
            [h2_vect, ~, ~] = easy_state(high2_mat,diag(B),eye(nNode),zeros(nNode));
            high2_states(:,i,k) = h2_vect;
            [h3_vect, ~, ~] = easy_state(high3_mat,diag(B),eye(nNode),zeros(nNode));
            high3_states(:,i,k) = h3_vect;
            [l_vect, ~, ~, l_v2] = easy_state(low_mat,diag(B),eye(nNode),zeros(nNode));
            low_states(:,i,k) = l_vect;
            low_states_ev2(:,i,k) = l_v2;
            
            z_vect_all = zeros(nNode, nZero);
            z_ev2_all = zeros(nNode, nZero);
            for j = 1:nZero
                [z_vect, ~, ~, z_ev2] = easy_state(zero_mat(:,:,j),diag(B),eye(nNode),zeros(nNode));
                z_vect_all(:,j) = z_vect;
                z_ev2_all(:,j) = z_ev2;
            end
            zero_states(:,i,k) = mean(z_vect_all,2);
            zero_states_ev2(:,i,k) = mean(z_ev2_all,2);
        end
    end
end

%% visualize


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
    
    plot_data.powspctrm = mean(high2_states(:,:,i),2);
    figure(3); clf
    ft_topoplotER(cfg,plot_data); colorbar; caxis([0.06,.13])
    saveas(gca, [img_dir, bands{i}, '_avg_high2_state_pr.png'], 'png')
    
    plot_data.powspctrm = mean(high3_states(:,:,i),2);
    figure(4); clf
    ft_topoplotER(cfg,plot_data); colorbar; caxis([0.06,.13])
    saveas(gca, [img_dir, bands{i}, '_avg_high3_state_pr.png'], 'png')
    
    plot_data.powspctrm = mean(high_states_ev2(:,:,i),2);
    figure(5); clf
    ft_topoplotER(cfg,plot_data); colorbar; caxis([-.1,.1])
    saveas(gca, [img_dir, bands{i}, '_avg_high_ev2_state_pr.png'], 'png')
    
    plot_data.powspctrm = mean(low_states_ev2(:,:,i),2);
    figure(6); clf
    ft_topoplotER(cfg,plot_data); colorbar; caxis([-.1,.1])
    saveas(gca, [img_dir, bands{i}, '_avg_low_ev2_state_pr.png'], 'png')
    
        plot_data.powspctrm = mean(zero_states(:,:,i),2);
    figure(7); clf
    ft_topoplotER(cfg,plot_data); colorbar; caxis([0.06,.13])
    saveas(gca, [img_dir, bands{i}, '_avg_zero_state.png'], 'png')
    
    plot_data.powspctrm = mean(zero_states_ev2(:,:,i),2);
    figure(8); clf
    ft_topoplotER(cfg,plot_data); colorbar; caxis([-.1,.1])
    saveas(gca, [img_dir, bands{i}, '_avg_zero_ev2_state.png'], 'png')
end

%% Categorize by lobe


h_region = [];
l_region = [];
l_region2 = [];
h_region2 = [];
h2_region = [];
h3_region = [];
z_region = [];
z_region2 = [];
region_ord = {};
band_ord = {};

cnt = 1;
for i = 1:numel(regions)
    for j = 1:numel(bands)
        load([top_dir, 'montages/', regions{i}, '_idx.mat'])
        idx = logical(idx);
        % get the mean edge between a lobes (total edges/number of edges), divided by
        % the total edges
        
        % stupid frmatting thing for R
        for k = 0:19
            region_ord{cnt+k} = regions{i};
            band_ord{cnt+k} = bands{j};
        end
        h_region = [h_region, mean(high_states(idx,:,j))];
        l_region = [l_region, mean(low_states(idx,:,j))];
        h2_region = [h2_region, mean(high2_states(idx,:,j))];
        h3_region = [h3_region, mean(high3_states(idx,:,j))];
        h_region2 = [h_region2, mean(high_states_ev2(idx,:,j))];
        l_region2 = [l_region2, mean(low_states_ev2(idx,:,j))];
        z_region = [z_region, mean(zero_states(idx,:,j))];
        z_region2 = [z_region2, mean(zero_states_ev2(idx,:,j))];
        
        cnt = cnt + nSubj;
    end
end

save([R_dir, 'grad/control_state_pr.mat'], 'region_ord', 'band_ord', 'h_region', 'l_region', 'h2_region', 'h3_region', 'h_region2', 'l_region2', 'z_region', 'z_region2')