%% Get some graph theory metrics for comparison
% look at path length, integration and segregation


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
nNode = 102;
nEdges = (nNode^2-nNode)/2;
load([save_dir, 'noise_sg.mat']);

% initialize
band_order = {};
subj_order = {};
high_eff = [];
high2_eff = [];
high3_eff = [];
low_eff = [];
zero_eff = [];
high_eff_rm = [];
high2_eff_rm = [];
high3_eff_rm = [];
low_eff_rm = [];
zero_eff_rm = [];
high_nc = [];
high2_nc = [];
high3_nc = [];
low_nc = [];
zero_nc = [];
high_pc = [];
high2_pc = [];
high3_pc = [];
low_pc = [];
zero_pc = [];

%% Get Motor index

load([top_dir, 'montages/Left_motor_idx.mat'])
lm = logical(idx);

% right motor seems like a stricter control
load([top_dir, 'montages/Right_motor_idx.mat'])
rm = logical(idx);

%% Loop through data

sens = sensors{1};
R_dir_s = [R_dir, sens, '/'];
cnt = 0; % for labels

for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        cnt = cnt + 1;
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
        
        %get matrix
        % scale to be stable
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
            zero_mat(:,:,j) = get_sg_matrix(nNode, subset(nbSG,:));
        end
        
        % get path length
        E_high = efficiency_wei(high_mat,2); % 2 for local efficiency
        E_high2 = efficiency_wei(high2_mat,2); % local efficiency
        E_high3 = efficiency_wei(high3_mat,2); % local efficiency
        E_low = efficiency_wei(low_mat,2); % local efficiency
        E_zero = zeros(nNode,nZero);
        for j = 1:nZero
            E_zero(:,j) = efficiency_wei(zero_mat(:,:,j),2);
        end
        
        % get only left motor
        E_high_lm = mean(E_high(lm));
        E_high2_lm = mean(E_high2(lm));
        E_high3_lm = mean(E_high3(lm));
        E_low_lm = mean(E_low(lm));
        E_zero_lm = mean(mean(E_zero(lm,:),1));
        
        % add to data structures
        band_order{cnt} = f;
        subj_order{cnt} = subj;
        high_eff = [high_eff; E_high_lm];
        high2_eff = [high2_eff; E_high2_lm];
        high3_eff = [high3_eff; E_high3_lm];
        low_eff = [low_eff; E_low_lm];
        zero_eff = [zero_eff; E_zero_lm];
        
        
        % repeat for rm
        E_high_rm = mean(E_high(rm));
        E_high2_rm = mean(E_high2(rm));
        E_high3_rm = mean(E_high3(rm));
        E_low_rm = mean(E_low(rm));
        E_zero_rm = mean(mean(E_zero(rm,:),1));
        
        % add to data structures
        band_order{cnt} = f;
        subj_order{cnt} = subj;
        high_eff_rm = [high_eff_rm; E_high_rm];
        high2_eff_rm = [high2_eff_rm; E_high2_rm];
        high3_eff_rm = [high3_eff_rm; E_high3_rm];
        low_eff_rm = [low_eff_rm; E_low_rm];
        zero_eff_rm = [zero_eff_rm; E_zero_rm];
        
        
        % now do some community detection
        M_high = community_louvain(high_mat,1); % gamma = 1
        M_high2 = community_louvain(high2_mat,1);
        M_high3 = community_louvain(high3_mat,1); 
        M_low = community_louvain(low_mat,1); 
        M_zero = zeros(nNode,nZero);
        for j = 1:nZero
            M_zero(:,j) = community_louvain(zero_mat(:,:,j),1);
        end
        
        % get the number of communities
        high_nc = [high_nc; numel(unique(M_high(lm)))];
        high2_nc = [high2_nc; numel(unique(M_high2(lm)))];
        high3_nc = [high3_nc; numel(unique(M_high3(lm)))];
        low_nc = [low_nc; numel(unique(M_low(lm)))];
        tmp = zeros(1,nZero);
        for j = 1:nZero
            tmp(j) = numel(unique(M_zero(:,j)));
        end
        zero_nc = [zero_nc; mean(tmp)];
        
        % participation coeff
        pc_high = participation_coef(high_mat,M_high); % gamma = 1
        pc_high2 = participation_coef(high2_mat,M_high2);
        pc_high3 = participation_coef(high3_mat,M_high3); 
        pc_low = participation_coef(low_mat,M_low); 
        pc_zero = zeros(nNode,nZero);
        for j = 1:nZero
            pc_zero(:,j) = participation_coef(zero_mat(:,:,j),M_zero(:,j));
        end
        
        high_pc = [high_pc; mean(pc_high(lm))];
        high2_pc = [high2_pc; mean(pc_high2(lm))];
        high3_pc = [high3_pc; mean(pc_high3(lm))];
        low_pc = [low_pc; mean(pc_low(lm))];
        zero_pc = [zero_pc; mean(mean(pc_zero(lm)))];
        
    end
end

% save
save([R_dir, 'grad/graph_stats.mat'], 'band_order', 'subj_order', 'high_eff', 'high2_eff', 'high3_eff', 'low_eff', 'zero_eff', 'high_eff_rm',...
    'high2_eff_rm', 'high3_eff_rm', 'low_eff_rm', 'zero_eff_rm', 'high_nc', 'high2_nc', 'high3_nc', 'low_nc', 'zero_nc', 'high_pc', 'high2_pc',...
    'high3_pc', 'low_pc', 'zero_pc')

