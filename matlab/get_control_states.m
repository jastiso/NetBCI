%% Get largest eigenvalue state
% get the "easiest to get to" state from the max and min
% also get the state associated with the smallest eigenvalue
% also get the easiest to reach for the second largest behavioral loading
% and get the state associated with the second largest eigenvalue
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
%regions = [{'Left_motor'}, {'Right_motor'}];


nNode = 102;
nEdges = (nNode^2-nNode)/2;
load([save_dir, 'noise_sg.mat']);

% initialize
high_states = zeros(nNode, nSubj, numel(bands));
highish_states = zeros(nNode, nSubj, numel(bands));
high_states_ev2 = zeros(nNode, nSubj, numel(bands));
low_states = zeros(nNode, nSubj, numel(bands));
small_states = zeros(nNode, nSubj, numel(bands));

%% Make control set

% get labels
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

load([top_dir, 'montages/Left_motor_idx.mat'])
B = idx;

% Get states for right temporal control state
nControl = sum(B);

% chosing right temporal
%load([top_dir, 'montages/Right_temporal_idx.mat'])
%B_control = idx;
%load([top_dir, 'montages/Right_occipital_idx.mat'])
% get two regions from occipital that are close
%B_control(97:98) = idx(97:98);
%B_control = double(B_control);
%sum(B_control)

% right motor seems like a stricter control
load([top_dir, 'montages/Right_motor_idx.mat'])
B_control = idx;

%% Loop through data

% this used to index over mag and grad, but we dont look at mag
sens = sensors{1};
R_dir_s = [R_dir, sens, '/'];

for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        f = bands{k};
        
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
        
        b_exp = subset(:,end);
        [~,bSG] = max(b_exp);
        [~, bSG2] = max(b_exp(b_exp<max(b_exp))); % second biggest element
        [~, nzSG] = min(nonzeros((b_exp))); % smallest nonzero - this is how we will opporationalize "low"
        
        % if you actualy want to look at 0s
        %if sum(b_exp == 0) <= 1
        %    [~,nbSG] = min(b_exp);
        %else
        %    nbSG = find(b_exp == 0);
        %end
        %nZero = numel(nbSG);
        
        
        
        
        %get gramian
        % scale to be stable
        high_mat = get_sg_matrix(nNode, subset(bSG,:));
        high_mat = high_mat./eigs(high_mat,1) - eye(nNode).*1.001;
        highish_mat = get_sg_matrix(nNode, subset(bSG2,:));
        highish_mat = highish_mat./eigs(highish_mat,1) - eye(nNode).*1.001;
        low_mat = get_sg_matrix(nNode, subset(nzSG,:));
        low_mat = low_mat./eigs(low_mat,1) - eye(nNode).*1.001;
        %low_mat = zeros(nNode, nNode, nZero);
        %for j = 1:nZero
        %    low_mat(:,:,j) = get_sg_matrix(nNode, subset(nbSG,:));
        %    low_mat(:,:,j) = low_mat(:,:,j)./eigs(low_mat(:,:,j),1) - eye(nNode).*1.001;
        %end
        
        
        % get largest eigenvector of grammian
        [h_vect, ~, small_vect, vect2] = easy_state(high_mat,diag(B),eye(nNode),zeros(nNode));
        high_states(:,i,k) = h_vect;
        small_states(:,i,k) = small_vect;
        high_states_ev2(:,i,k) = vect2;
        [hish_vect, ~, ~] = easy_state(highish_mat,diag(B),eye(nNode),zeros(nNode));
        highish_states(:,i,k) = hish_vect;
        [l_vect, ~, ~] = easy_state(low_mat,diag(B),eye(nNode),zeros(nNode));
        low_states(:,i,k) = l_vect;
        
        %l_vect_all = zeros(nNode, nZero);
        %for j = 1:nZero
        %    [l_vect, l_val] = easy_state(low_mat(:,:,j),diag(B),eye(nNode),zeros(nNode));
        %     l_vect_all(:,j) = l_vect;
        %end
        %low_states(:,i,k) = mean(l_vect_all,2);
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
    ft_topoplotER(cfg,plot_data); colorbar; caxis([0.05,.15])
    saveas(gca, [img_dir, bands{i}, '_avg_high_state.png'], 'png')
    
    plot_data.powspctrm = mean(low_states(:,:,i),2);
    figure(2); clf
    ft_topoplotER(cfg,plot_data); colorbar; caxis([0.05,.15])
    saveas(gca, [img_dir, bands{i}, '_avg_low_state.png'], 'png')
    
    plot_data.powspctrm = mean(small_states(:,:,i),2);
    figure(3); clf
    ft_topoplotER(cfg,plot_data); colorbar; 
    saveas(gca, [img_dir, bands{i}, '_avg_small_state.png'], 'png')
    
    plot_data.powspctrm = mean(highish_states(:,:,i),2);
    figure(4); clf
    ft_topoplotER(cfg,plot_data); colorbar; caxis([0.05,.15])
    saveas(gca, [img_dir, bands{i}, '_avg_highish_state.png'], 'png')
    
    plot_data.powspctrm = mean(high_states_ev2(:,:,i),2);
    figure(5); clf
    ft_topoplotER(cfg,plot_data); colorbar; 
    saveas(gca, [img_dir, bands{i}, '_avg_ev2_state.png'], 'png')
    
end

%% Categorize by lobe


h_region = [];
l_region = [];
h_region2 = [];
hish_region = [];
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
        hish_region = [hish_region, mean(highish_states(idx,:,j))];
        h_region2 = [h_region2, mean(high_states_ev2(idx,:,:))];
        cnt = cnt + nSubj;
    end
end

save([R_dir, 'grad/control_state.mat'], 'region_ord', 'band_ord', 'h_region', 'l_region', 'hish_region', 'h_region2')



%% Other control set - Loop through data --------------------------------

cont_high_states = zeros(nNode, nSubj, numel(bands));
cont_high_states_ev2 = zeros(nNode, nSubj, numel(bands));
cont_low_states = zeros(nNode, nSubj, numel(bands));
cont_low_states_ev2 = zeros(nNode, nSubj, numel(bands));

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
            [cont_h_vect, h_val, ~, cont_h_vect2] = easy_state(high_mat,diag(B_control),eye(nNode),zeros(nNode));
            [cont_l_vect, l_val, ~, cont_l_vect2] = easy_state(low_mat,diag(B_control),eye(nNode),zeros(nNode));
            
            cont_high_states(:,i,k) = cont_h_vect;
            cont_low_states(:,i,k) = cont_l_vect;
            cont_high_states_ev2(:,i,k) = cont_h_vect2;
            cont_low_states_ev2(:,i,k) = cont_l_vect2;
        end
    end
end

%% visualize


cfg.layout = [top_dir, 'layouts/neuromag306cmb.lay'];
fid = fopen([top_dir, 'layouts/neuromag306cmb.txt'], 'r');

for i = 1:numel(bands)
    plot_data.powspctrm = mean(cont_high_states(:,:,i),2);
    plot_data.label = cmb_labels;
    plot_data.dimord = 'chan_freq';
    plot_data.freq = mean(freqs(i,:));
    plot_data.cfg = [];
    figure(1); clf
    ft_topoplotER(cfg,plot_data); colorbar
    saveas(gca, [img_dir, bands{i}, '_cont_avg_high_state.png'], 'png')
    
    plot_data.powspctrm = mean(cont_low_states(:,:,i),2);
    figure(1); clf
    ft_topoplotER(cfg,plot_data); colorbar
    saveas(gca, [img_dir, bands{i}, '_cont_avg_low_state.png'], 'png')
    
    % high
    plot_data.powspctrm = mean(cont_high_states(:,:,i),2);
    figure(2); clf
    ft_topoplotER(cfg,plot_data); colorbar
    saveas(gca, [img_dir, bands{i}, '_cont_avg_high_state.png'], 'png')
    
    % ev2
    plot_data.powspctrm = mean(cont_high_states_ev2(:,:,i),2);
    figure(3); colorbar
    ft_topoplotER(cfg,plot_data); colorbar
    saveas(gca, [img_dir, bands{i}, '_cont_high_ev2_state.png'], 'png')
    
     % ev2
    plot_data.powspctrm = mean(cont_low_states_ev2(:,:,i),2);
    figure(4); colorbar
    ft_topoplotER(cfg,plot_data); colorbar
    saveas(gca, [img_dir, bands{i}, '_cont_low_ev2_state.png'], 'png')
end

%% Categorize by lobe


cont_h_region = [];
cont_l_region = [];
cont_h_region2 = [];
cont_l_region2 = [];
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
        cont_h_region = [cont_h_region, mean(cont_high_states(idx,:,j))];
        cont_l_region = [cont_l_region, mean(cont_low_states(idx,:,j))];
        cont_h_region2 = [cont_h_region2, mean(cont_high_states_ev2(idx,:,j))];
        cont_l_region2 = [cont_l_region2, mean(cont_low_states_ev2(idx,:,j))];
        
        cnt = cnt + nSubj;
    end
end

save([R_dir, 'grad/cont_control_state.mat'], 'region_ord', 'band_ord', 'cont_h_region', 'cont_l_region', 'cont_h_region2', 'cont_l_region2')