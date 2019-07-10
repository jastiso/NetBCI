%% Consistency thresholds
% look at how consistency shanges for different thresholds
clear all

%% Global Variables

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath(genpath('/Users/stiso/Documents/MATLAB/BCT/'))

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

    end
end


%% make and plot consensus - consistency

thr = .1:.1:.9;
sum_high = zeros(numel(thr), numel(bands));
sum_high2 = zeros(numel(thr), numel(bands));
sum_high3 = zeros(numel(thr), numel(bands));
sum_low = zeros(numel(thr), numel(bands));
sp_high = zeros(numel(thr), numel(bands));
sp_high2 = zeros(numel(thr), numel(bands));
sp_high3 = zeros(numel(thr), numel(bands));
sp_low = zeros(numel(thr), numel(bands));

for i = 1:numel(bands)
    curr_high = consensus_high(:,:,:,i);
    curr_high2 = consensus_high2(:,:,:,i);
    curr_high3 = consensus_high3(:,:,:,i);
    curr_low = consensus_low(:,:,:,i);
    
    for j = 1:numel(thr)
        c_high = zeros(nNode,nNode, nSubj);
        c_high2 = zeros(nNode,nNode, nSubj);
        c_high3 = zeros(nNode,nNode, nSubj);
        c_low = zeros(nNode,nNode, nSubj);

        % threshold
        for k = 1:nSubj
            c_high(:,:,k) = thresh_mat(curr_high(:,:,k),thr(j));
            c_high2(:,:,k) = thresh_mat(curr_high2(:,:,k),thr(j));
            c_high3(:,:,k) = thresh_mat(curr_high3(:,:,k),thr(j));
            c_low(:,:,k) = thresh_mat(curr_low(:,:,k),thr(j));
        end
        
        % get consistency
        c_high_avg = sum(c_high ~=0, 3); % ~=0 for consistency
        c_high2_avg = sum(c_high2 ~=0, 3);
        c_high3_avg = sum(c_high3 ~=0, 3);
        c_low_avg = sum(c_low ~=0, 3);
        
        % get sparsity
        sp_high(j,i) = mean(squeeze(sum(sum(c_high == 0))./numel(c_high_avg))); % avg sparsity across all subjs
        sp_high2(j,i) = mean(squeeze(sum(sum(c_high2 == 0))./numel(c_high2_avg))); % avg sparsity across all subjs
        sp_high3(j,i) = mean(squeeze(sum(sum(c_high3 == 0))./numel(c_high3_avg))); % avg sparsity across all subjs
        sp_low(j,i) = mean(squeeze(sum(sum(c_low == 0))./numel(c_low_avg))); % avg sparsity across all subjs
        
        % get number of edges remaining and max
        sum_high(j,i) = sum(sum(c_high_avg))/numel(c_high_avg); % same number in each sg
        sum_high2(j,i) = sum(sum(c_high2_avg))/numel(c_high_avg); % same number in each sg
        sum_high3(j,i) = sum(sum(c_high3_avg))/numel(c_high_avg); % same number in each sg
        sum_low(j,i) = sum(sum(c_low_avg))/numel(c_high_avg); % same number in each sg
       
        
    end
end

%% Null data

load([save_dir, 'pr_noise_sg.mat']);

sens = sensors{1};

for i = subjs
    subj = sprintf('%03d', i);
    load([top_dir, 'idx_', sens, '.mat'])
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nSubj %03d... band %s...', i, f)
        
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

    end
end

sum_high_pr = zeros(numel(thr), numel(bands));
sum_high2_pr = zeros(numel(thr), numel(bands));
sum_high3_pr = zeros(numel(thr), numel(bands));
sum_low_pr = zeros(numel(thr), numel(bands));

sp_high_pr = zeros(numel(thr), numel(bands));
sp_high2_pr = zeros(numel(thr), numel(bands));
sp_high3_pr = zeros(numel(thr), numel(bands));
sp_low_pr = zeros(numel(thr), numel(bands));

for i = 1:numel(bands)
    curr_high = consensus_high(:,:,:,i);
    curr_high2 = consensus_high2(:,:,:,i);
    curr_high3 = consensus_high3(:,:,:,i);
    curr_low = consensus_low(:,:,:,i);
    
    for j = 1:numel(thr)
        c_high = zeros(nNode,nNode, nSubj);
        c_high2 = zeros(nNode,nNode, nSubj);
        c_high3 = zeros(nNode,nNode, nSubj);
        c_low = zeros(nNode,nNode, nSubj);

        % threshold
        for k = 1:nSubj
            c_high(:,:,k) = thresh_mat(curr_high(:,:,k),thr(j));
            c_high2(:,:,k) = thresh_mat(curr_high2(:,:,k),thr(j));
            c_high3(:,:,k) = thresh_mat(curr_high3(:,:,k),thr(j));
            c_low(:,:,k) = thresh_mat(curr_low(:,:,k),thr(j));
        end
        % get consistency
        c_high_avg = sum(c_high ~=0, 3); % ~=0 for consistency
        c_high2_avg = sum(c_high2 ~=0, 3);
        c_high3_avg = sum(c_high3 ~=0, 3);
        c_low_avg = sum(c_low ~=0, 3);        
        
        % get sparsity
        sp_high_pr(j,i) = mean(squeeze(sum(sum(c_high == 0))./numel(c_high_avg))); % avg sparsity across all subjs
        sp_high2_pr(j,i) = mean(squeeze(sum(sum(c_high2 == 0))./numel(c_high2_avg))); % avg sparsity across all subjs
        sp_high_pr(j,i) = mean(squeeze(sum(sum(c_high3 == 0))./numel(c_high3_avg))); % avg sparsity across all subjs
        sp_low_pr(j,i) = mean(squeeze(sum(sum(c_low == 0))./numel(c_low_avg))); % avg sparsity across all subjs
        
        % get number of edges remaining
        sum_high_pr(j,i) = sum(sum(c_high_avg))/numel(c_high_avg); % same number in each sg
        sum_high2_pr(j,i) = sum(sum(c_high2_avg))/numel(c_high_avg); % same number in each sg
        sum_high3_pr(j,i) = sum(sum(c_high3_avg))/numel(c_high_avg); % same number in each sg
        sum_low_pr(j,i) = sum(sum(c_low_avg))/numel(c_high_avg); % same number in each sg
        
   
    end
end

%% Plot consistency for all thresh

gamma = [232/255,200/255,134/255];
alpha = [0/255,89/255,115/255];
beta = [151/255, 0/255,0/255];

figure(1); clf
plot(sum_high(:,1), 'color', alpha, "linewidth", 3); hold on
plot(sum_high(:,2), 'color', beta, "linewidth", 3)
plot(sum_high(:,3), 'color', gamma, "linewidth", 3)
plot(sum_high_pr(:,1), 'color', alpha, "linewidth", 3, 'linestyle', '--')
plot(sum_high_pr(:,2), 'color', beta, "linewidth", 3, 'linestyle', '--')
plot(sum_high_pr(:,3), 'color', gamma, "linewidth", 3, 'linestyle', '--')
title('high')
saveas(gca, [img_dir, 'cinstistency_thresh_high.png'], 'png')

figure(2); clf
plot(sum_high2(:,1), 'color', alpha, "linewidth", 3); hold on
plot(sum_high2(:,2), 'color', beta, "linewidth", 3)
plot(sum_high2(:,3), 'color', gamma, "linewidth", 3)
plot(sum_high2_pr(:,1), 'color', alpha, "linewidth", 3, 'linestyle', '--')
plot(sum_high2_pr(:,2), 'color', beta, "linewidth", 3, 'linestyle', '--')
plot(sum_high2_pr(:,3), 'color', gamma, "linewidth", 3, 'linestyle', '--')
title('high2')
saveas(gca, [img_dir, 'cinstistency_thresh_high2.png'], 'png')

figure(3); clf
plot(sum_high3(:,1), 'color', alpha, "linewidth", 3); hold on
plot(sum_high3(:,2), 'color', beta, "linewidth", 3)
plot(sum_high3(:,3), 'color', gamma, "linewidth", 3)
plot(sum_high3_pr(:,1), 'color', alpha, "linewidth", 3, 'linestyle', '--')
plot(sum_high3_pr(:,2), 'color', beta, "linewidth", 3, 'linestyle', '--')
plot(sum_high3_pr(:,3), 'color', gamma, "linewidth", 3, 'linestyle', '--')
title('high3')
saveas(gca, [img_dir, 'cinstistency_thresh_high3.png'], 'png')

figure(4); clf
plot(sum_low(:,1), 'color', alpha, "linewidth", 3); hold on
plot(sum_low(:,2), 'color', beta, "linewidth", 3)
plot(sum_low(:,3), 'color', gamma, "linewidth", 3)
plot(sum_low_pr(:,1), 'color', alpha, "linewidth", 3, 'linestyle', '--')
plot(sum_low_pr(:,2), 'color', beta, "linewidth", 3, 'linestyle', '--')
plot(sum_low_pr(:,3), 'color', gamma, "linewidth", 3, 'linestyle', '--')
title('low')
saveas(gca, [img_dir, 'cinstistency_thresh_low.png'], 'png')


%% Plto sparsity


gamma = [232/255,200/255,134/255];
alpha = [0/255,89/255,115/255];
beta = [151/255, 0/255,0/255];

figure(1); clf
plot(sum_high(:,1) - sum_high_pr(:,1), 'color', alpha, "linewidth", 3); hold on
plot(sum_high(:,2) - sum_high_pr(:,2), 'color', beta, "linewidth", 3)
plot(sum_high(:,3) - sum_high_pr(:,3), 'color', gamma, "linewidth", 3)
% plot(sp_high_pr(:,1), 'color', alpha, "linewidth", 3, 'linestyle', '--')
% plot(sp_high_pr(:,2), 'color', beta, "linewidth", 3, 'linestyle', '--')
% plot(sp_high_pr(:,3), 'color', gamma, "linewidth", 3, 'linestyle', '--')
title('high')
saveas(gca, [img_dir, 'cinstistency_sp_high.png'], 'png')

figure(2); clf
plot(sum_high2(:,1) - sum_high2_pr(:,1), 'color', alpha, "linewidth", 3); hold on
plot(sum_high2(:,2) - sum_high2_pr(:,2), 'color', beta, "linewidth", 3)
plot(sum_high2(:,3) - sum_high2_pr(:,3), 'color', gamma, "linewidth", 3)
% plot(sp_high2_pr(:,1), 'color', alpha, "linewidth", 3, 'linestyle', '--')
% plot(sp_high2_pr(:,2), 'color', beta, "linewidth", 3, 'linestyle', '--')
% plot(sp_high2_pr(:,3), 'color', gamma, "linewidth", 3, 'linestyle', '--')
title('high2')
saveas(gca, [img_dir, 'cinstistency_sp_high2.png'], 'png')

figure(3); clf
plot(sum_high3(:,1) - sum_high3_pr(:,1), 'color', alpha, "linewidth", 3); hold on
plot(sum_high3(:,2) - sum_high3_pr(:,2), 'color', beta, "linewidth", 3)
plot(sum_high3(:,3) - sum_high3_pr(:,3), 'color', gamma, "linewidth", 3)
% plot(sp_high3_pr(:,1), 'color', alpha, "linewidth", 3, 'linestyle', '--')
% plot(sp_high3_pr(:,2), 'color', beta, "linewidth", 3, 'linestyle', '--')
% plot(sp_high3_pr(:,3), 'color', gamma, "linewidth", 3, 'linestyle', '--')
title('high3')
saveas(gca, [img_dir, 'cinstistency_sp_high3.png'], 'png')

figure(4); clf
plot(sum_low(:,1)  - sum_low_pr(:,1), 'color', alpha, "linewidth", 3); hold on
plot(sum_low(:,2) - sum_low_pr(:,2), 'color', beta, "linewidth", 3)
plot(sum_low(:,3) - sum_low_pr(:,3), 'color', gamma, "linewidth", 3)
%plot(sp_low_pr(:,2), 'color', beta, "linewidth", 3, 'linestyle', '--')
%plot(sp_low_pr(:,3), 'color', gamma, "linewidth", 3, 'linestyle', '--')
title('low')
saveas(gca, [img_dir, 'cinstistency_sp_low.png'], 'png')