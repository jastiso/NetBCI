%% Consensus Subgraphs
% make consesnsus subgraphs for each


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
load([save_dir, 'pr_noise_sg.mat']);

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
        %consensus_zero(:,:,i,k) = zero_mat;
    end
end


%% make and plot consensus - consistency

thr = [.75];
sum_high = zeros(numel(thr), numel(bands));
sum_high2 = zeros(numel(thr), numel(bands));
sum_high3 = zeros(numel(thr), numel(bands));
sum_low = zeros(numel(thr), numel(bands));

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
        c_high_avg = sum(c_high ~=0, 3); % ~=0 for consistency
        c_high2_avg = sum(c_high2 ~=0, 3);
        c_high3_avg = sum(c_high3 ~=0, 3);
        c_low_avg = sum(c_low ~=0, 3);
        %c_zero_avg = sum(c_zero ~= 0, 3);
        
        % get number of edges remaining
        sum_high(j,i) = sum(sum(c_high_avg))/numel(c_high_avg); % same number in each sg
        sum_high2(j,i) = sum(sum(c_high2_avg))/numel(c_high_avg); % same number in each sg
        sum_high3(j,i) = sum(sum(c_high3_avg))/numel(c_high_avg); % same number in each sg
        sum_low(j,i) = sum(sum(c_low_avg))/numel(c_high_avg); % same number in each sg
        
        
        % save consistent mats for strictest thresh
        save([save_dir, 'consistent_mats_', bands{i}, '_pr.mat'], 'c_high_avg', 'c_high2_avg', 'c_high3_avg', 'c_low_avg')
        
        % add some additional thresholding
        c_high_avg(c_high_avg <= 5) = 0;
        c_high2_avg(c_high2_avg <= 5) = 0;
        c_high3_avg(c_high3_avg <= 5) = 0;
        c_low_avg(c_low_avg <= 5) = 0;
        
        % save stricter thresholding to report later
        [r,c] = find(c_high_avg > 9);
        c_high_max{i,j} = [r,c,c_high_avg(r,c)];
        [r,c] = find(c_high2_avg > 9);
        c_high2_max{i,j} = [r,c,c_high2_avg(r,c)];
        [r,c] = find(c_high3_avg > 9);
        c_high3_max{i,j} = [r,c,c_high3_avg(r,c)];
        [r,c] = find(c_low_avg > 9);
        c_low_max{i,j} = [r,c,c_low_avg(r,c)];
        
        % cluster
        [ci_high,q_high] = community_louvain(c_high_avg,1);
        [~,idx_high] = sort(ci_high);
        [ci_high2,q_high2] = community_louvain(c_high2_avg,1);
        [~,idx_high2] = sort(ci_high2);
        [ci_high3,q_high3] = community_louvain(c_high3_avg,1);
        [~,idx_high3] = sort(ci_high3);
        [ci_low,q_low] = community_louvain(c_low_avg,1);
        [~,idx_low] = sort(ci_low);
        
        
        % save for gephi - coordinates, and community in CSV, and matrix in
        % txt
%         dlmwrite(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_high_matrix_pr.txt'],c_high_avg, 'delimiter', '\t')
%         dlmwrite(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_high2_matrix_pr.txt'],c_high2_avg, 'delimiter', '\t')
%         dlmwrite(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_high3_matrix_pr.txt'],c_high3_avg, 'delimiter', '\t')
%         dlmwrite(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_low_matrix_pr.txt'],c_low_avg, 'delimiter', '\t')
%         % csv
%         fid = fopen([top_dir, 'layouts/neuromag306cmb.txt'], 'r');
%         coord = textscan(fid, '%d %10f %10f %d %d %s', 'Delimiter', '\n');
%         fclose(fid);
%         header = ['x,y,z,community,'];
%         
%         fid = fopen(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_high_attributes_pr.csv'],'w');
%         fprintf(fid,'%s\n',header); fclose(fid);
%         dlmwrite(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_high_attributes_pr.csv'],[coord{2}, coord{3}, ones(nNode,1), ci_high],'-append');
%         
%         fid = fopen(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_high2_attributes_pr.csv'],'w');
%         fprintf(fid,'%s\n',header); fclose(fid);
%         dlmwrite(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_high2_attributes_pr.csv'],[coord{2}, coord{3}, ones(nNode,1), ci_high2],'-append');
%         
%         fid = fopen(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_high3_attributes_pr.csv'],'w');
%         fprintf(fid,'%s\n',header); fclose(fid);
%         dlmwrite(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_high3_attributes_pr.csv'],[coord{2}, coord{3}, ones(nNode,1), ci_high3],'-append');
%         
%         fid = fopen(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_low_attributes_pr.csv'],'w');
%         fprintf(fid,'%s\n',header); fclose(fid);
%         dlmwrite(['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/wpli/gephi/', bands{i}, '_', num2str(j), '_low_attributes_pr.csv'],[coord{2}, coord{3}, ones(nNode,1), ci_low],'-append');
        
    end
end

save([save_dir, 'consistent_sums_pr.mat'], 'sum_high', 'sum_high2', 'sum_high3', 'sum_low')
save([save_dir, 'most_consistent_pr.mat'], 'c_high_max', 'c_high2_max', 'c_high3_max', 'c_low_max')