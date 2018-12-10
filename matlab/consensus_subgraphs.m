%% Consensus Subgraphs
% make consesnsus subgraphs for each


%% Global Variables

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/param/';
save_dir = [top_dir, 'GroupAvg/wpli/analysis/'];
img_dir = [top_dir, 'GroupAvg/wpli/images/'];

subjs = [1:20];
nSubj = numel(subjs);
nNode = 102;
nEdge = (nNode^2-nNode)/2;
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'}];
sensors = [{'grad'}];
load([save_dir, 'noise_sg.mat']);

% initialize
consensus_high = zeros(nNode,nNode, nSubj, numel(bands));
consensus_high2 = zeros(nNode,nNode, nSubj, numel(bands));
consensus_high3 = zeros(nNode,nNode, nSubj, numel(bands));
consensus_low = zeros(nNode,nNode, nSubj, numel(bands));
consensus_zero = zeros(nNode,nNode, nSubj, numel(bands));

% get labels for topoplot
labels = [];
load(['/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/Session1/1_test01/RawData_MEG_Subj001_Ses1_test01.mat'])
data = RawData_MEG_Subj001_Ses1_test01;
clear RawData_MEG_Subj001_Ses1_test01

sens_idx = strcmp(data.meg_sensors.chantype,'megplanar');
grad_label = data.label(sens_idx);
cnt = 1;
for n = 1:2:sum(sens_idx)
    labels{cnt,1} = [grad_label{n}, '+', grad_label{n+1}(end-3:end)];
    cnt = cnt + 1;
end

%% Loop through data

sens = sensors{1};

for i = subjs
    subj = sprintf('%03d', i);
    load([top_dir, 'idx_', sens, '.mat'])
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nSubj %03d... band %s...', i, f)
        
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
        [~,tmp_idx] = sort(b_exp,'descend');
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
        
        % add to consensus
        % 1st
        high_mat = get_sg_matrix(nNode, subset(bSG,:));
        %2nd
        high2_mat = get_sg_matrix(nNode, subset(bSG2,:));
        %3rd
        high3_mat = get_sg_matrix(nNode, subset(bSG3,:));
        %lowest
        low_mat = get_sg_matrix(nNode, subset(nzSG,:));
        zero_mat = zeros(nNode, nNode, nZero);
        for j = 1:nZero
            zero_mat(:,:,j) = get_sg_matrix(nNode, subset(nbSG(j),:));
        end
        
        consensus_high(:,:,i,k) = high_mat;
        consensus_high2(:,:,i,k) = high2_mat;
        consensus_high3(:,:,i,k) = high3_mat;
        consensus_low(:,:,i,k) = low_mat;
        %consensus_zero(:,:,i,k) = zero_mat;
    end
end


%% make and plot consensus - consistency

thr = [.25 .5 .75];

for i = 1:numel(bands)
    curr_high = consensus_high(:,:,:,i);
    curr_high2 = consensus_high2(:,:,:,i);
    curr_high3 = consensus_high3(:,:,:,i);
    curr_low = consensus_low(:,:,:,i);
    %curr_zero = consensus_zero(:,:,:,i);
    
    for j = 1:numel(thr)
        c_high = zeros(nNode,nNode, nSubj);
        c_high2 = zeros(nNode,nNode, nSubj);
        c_high3 = zeros(nNode,nNode, nSubj);
        c_low = zeros(nNode,nNode, nSubj);
        %c_zero = zeros(nNode,nNode, nSubj);
        % threshold
        for k = 1:nSubj
            c_high(:,:,k) = thresh_mat(curr_high(:,:,k),thr(j));
            c_high2(:,:,k) = thresh_mat(curr_high2(:,:,k),thr(j));
            c_high3(:,:,k) = thresh_mat(curr_high3(:,:,k),thr(j));
            c_low(:,:,k) = thresh_mat(curr_low(:,:,k),thr(j));
            %c_zero(:,:,k) = thresh_mat(curr_zero(:,:,k),thr(j));
        end
        % get consistency
        c_high_avg = sum(c_high, 3); % ~=0 for consistency
        c_high2_avg = sum(c_high2, 3);
        c_high3_avg = sum(c_high3, 3);
        c_low_avg = sum(c_low, 3);
        %c_zero_avg = sum(c_zero ~= 0, 3);
        
        % cluster
        [ci_high,q_high] = community_louvain(c_high_avg,1);
        [~,idx_high] = sort(ci_high);
        [ci_high2,q_high2] = community_louvain(c_high2_avg,1);
        [~,idx_high2] = sort(ci_high2);
        [ci_high3,q_high3] = community_louvain(c_high3_avg,1);
        [~,idx_high3] = sort(ci_high3);
        [ci_low,q_low] = community_louvain(c_low_avg,1);
        [~,idx_low] = sort(ci_low);
        
        %plot
        figure(1); clf
        imagesc(c_high_avg(idx_high,idx_high)); colorbar; colormap('pink') % change the colormap
        %caxis([0,20])
        %saveas(gca,[img_dir, 'consesnus_high_', bands{i}, '_', num2str(j), '.png'], 'png');
        
        figure(2); clf
        imagesc(c_high2_avg(idx_high2,idx_high2)); colorbar; colormap('pink')
        %caxis([0,20])
        %saveas(gca,[img_dir, 'consesnus_high2_', bands{i}, '_', num2str(j), '.png'], 'png');
        
        figure(3); clf
        imagesc(c_high3_avg(idx_high3,idx_high3)); colorbar; colormap('pink')
        %caxis([0,20])
        %saveas(gca,[img_dir, 'consesnus_high3_', bands{i}, '_', num2str(j), '.png'], 'png');
        
        figure(4); clf
        imagesc(c_low_avg(idx_low,idx_low)); colorbar; colormap('pink')
        %caxis([0,20])
        %saveas(gca,[img_dir, 'consesnus_low_', bands{i}, '_', num2str(j), '.png'], 'png');
        
        %figure(5); clf
        %imagesc(c_zero_avg); colorbar
        %caxis([0,20])
        %saveas(gca,[img_dir, 'consesnus_zero_', bands{i}, '_', num2str(j), '.png'], 'png');
        
        % get into fieldtrip format for topoplot
        cfg = [];
        cfg.style = 'straight';
        cfg.interpolatenans = 'no';
        %cfg.interplimits = 'electrodes';
        cfg.layout = [top_dir, 'layouts/neuromag306cmb.lay'];
        
        plot_data.powspctrm = ci_high;
        plot_data.label = labels;
        plot_data.dimord = 'chan_freq';
        plot_data.freq = 6;
        figure(1); clf
        ft_topoplotER(cfg,plot_data); colorbar
        % so far it was the same as above, now change the colormap
        ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
        colormap(flipud(brewermap(numel(unique(ci_high)),'Set2'))) % change the colormap
        saveas(gca, [img_dir, 'consesnus_high_', bands{i}, '_', num2str(j), '_topo.png'], 'png')
        
        plot_data.powspctrm = ci_high2;
        figure(2); clf
        ft_topoplotER(cfg,plot_data); colorbar
        % so far it was the same as above, now change the colormap
        ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
        colormap(flipud(brewermap(numel(unique(ci_high2)),'Set2'))) % change the colormap
        saveas(gca, [img_dir, 'consesnus_high2_', bands{i}, '_', num2str(j), '_topo.png'], 'png')
        
        plot_data.powspctrm = ci_high3;
        figure(3); clf
        ft_topoplotER(cfg,plot_data); colorbar
        % so far it was the same as above, now change the colormap
        ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
        colormap(flipud(brewermap(numel(unique(ci_high3)),'Set2'))) % change the colormap
        saveas(gca, [img_dir, 'consesnus_high3_', bands{i}, '_', num2str(j), '_topo.png'], 'png')
        
        plot_data.powspctrm = ci_low;
        figure(4); clf
        ft_topoplotER(cfg,plot_data); colorbar
        % so far it was the same as above, now change the colormap
        ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
        colormap(flipud(brewermap(numel(unique(ci_low)),'Set2'))) % change the colormap
        saveas(gca, [img_dir, 'consesnus_low_', bands{i}, '_', num2str(j), '_topo.png'], 'png')
        
    end
end
