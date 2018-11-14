%% Optimal control to motor imagery pattern

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath('/Users/stiso/Documents/MATLAB/control_example/')
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

% control parameters
rho = 1;
T = 1;
S = eye(nNode);
x0 = zeros(nNode,1);

% initialize
u_high = [];
u_high2 = [];
u_high3 = [];
u_low = [];
u_zero = [];
band_order = {};
subj_order = {};

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

% right motor seems like a stricter control
load([top_dir, 'montages/Right_motor_idx.mat'])
B_control = idx;

%% Make target states
% make sure these are in the order of bands

% alpha, suppression in left motor
xT(:,1) = -B;

% beta - contralateral suppression, and ipsilateral activation
xT(:,2) = -B + B_control;

% gamma...contralateral activation?
xT(:,3) = B;

%% Loop through data

% this used to index over mag and grad, but we dont look at mag
sens = sensors{1};
R_dir_s = [R_dir, sens, '/'];

cnt = 1; % for order
for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nSubj %03d... band %s...', i, f)

        % get labels
        band_order{cnt} = f;
        subj_order{cnt} = subj;
        cnt = cnt + 1;

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
        
        % optimal control
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

        % optim_fun(A, T, B, x0, xf, rho, S)
        [U_h, err_h] = get_opt_energy(high_mat, T, B, x0, xT(:,k), rho, S);
        u_high = [u_high; U_h];
        [U_h2, err_h2] = get_opt_energy(high2_mat, T, B, x0, xT(:,k), rho, S);
        u_high2 = [u_high2; U_h2];
        [U_h3, err_h3] = get_opt_energy(high3_mat, T, B, x0, xT(:,k), rho, S);
        u_high3 = [u_high3; U_h3];
        [U_l, err_l] = get_opt_energy(low_mat, T, B, x0, xT(:,k), rho, S);
        u_low = [u_low; U_l];
        U_z_all = zeros(nZero,1);
        for j = 1:nZero
            [U_z, err_z] = get_opt_energy(zero_mat(:,:,j), T, B, x0, xT(:,k), rho, S);
            U_z_all(j) = U_z;
        end
        u_zero = [u_zero; mean(U_z_all)];
    end
end

save([R_dir, 'grad/opt_energy.mat'], 'band_order', 'subj_order', 'u_high', 'u_high2', 'u_high3', 'u_low', 'u_zero')



%% Loop through data - other control set

% initialize
u_high_c = [];
u_high2_c = [];
u_high3_c = [];
u_low_c = [];
u_zero_c = [];
band_order_c = {};
subj_order_c = {};


% this used to index over mag and grad, but we dont look at mag
sens = sensors{1};
R_dir_s = [R_dir, sens, '/'];

cnt = 1; % for order
for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nControl Subj %03d... band %s...', i, f)

        % get labels
        band_order_c{cnt} = f;
        subj_order_c{cnt} = subj;
        cnt = cnt + 1;

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
        
        % optimal control
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

        % optim_fun(A, T, B, x0, xf, rho, S)
        [U_h, err_h] = get_opt_energy(high_mat, T, B_control, x0, xT(:,k), rho, S);
        u_high_c = [u_high_c; U_h];
        [U_h2, err_h2] = get_opt_energy(high2_mat, T, B_control, x0, xT(:,k), rho, S);
        u_high2_c = [u_high2_c; U_h2];
        [U_h3, err_h3] = get_opt_energy(high3_mat, T, B_control, x0, xT(:,k), rho, S);
        u_high3_c = [u_high3_c; U_h3];
        [U_l, err_l] = get_opt_energy(low_mat, T, B_control, x0, xT(:,k), rho, S);
        u_low_c = [u_low_c; U_l];
        U_z_all = zeros(nZero,1);
        for j = 1:nZero
            [U_z, err_z] = get_opt_energy(zero_mat(:,:,j), T, B_control, x0, xT(:,k), rho, S);
            U_z_all(j) = U_z;
        end
        u_zero_c = [u_zero_c; mean(U_z_all)];
    end
end

save([R_dir, 'grad/opt_energy_control.mat'], 'band_order_c', 'subj_order_c', 'u_high_c', 'u_high2_c', 'u_high3_c', 'u_low_c', 'u_zero_c')



%% Loop through data - PR

load([save_dir, 'pr_noise_sg.mat']);

% initialize
u_high_pr = [];
u_high2_pr = [];
u_high3_pr = [];
u_low_pr = [];
u_zero_pr = [];
band_order_pr = {};
subj_order_pr = {};


% this used to index over mag and grad, but we dont look at mag
sens = sensors{1};
R_dir_s = [R_dir, sens, '/'];

cnt = 1; % for order
for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nUPR Subj %03d... band %s...', i, f)

        % get labels
        band_order_pr{cnt} = f;
        subj_order_pr{cnt} = subj;
        cnt = cnt + 1;

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
        
        % optimal control
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

        % optim_fun(A, T, B, x0, xf, rho, S)
        [U_h, err_h] = get_opt_energy(high_mat, T, B, x0, xT(:,k), rho, S);
        u_high_pr = [u_high_pr; U_h];
        [U_h2, err_h2] = get_opt_energy(high2_mat, T, B, x0, xT(:,k), rho, S);
        u_high2_pr = [u_high2_pr; U_h2];
        [U_h3, err_h3] = get_opt_energy(high3_mat, T, B, x0, xT(:,k), rho, S);
        u_high3_pr = [u_high3_pr; U_h3];
        [U_l, err_l] = get_opt_energy(low_mat, T, B, x0, xT(:,k), rho, S);
        u_low_pr = [u_low_pr; U_l];
        U_z_all = zeros(nZero,1);
        for j = 1:nZero
            [U_z, err_z] = get_opt_energy(zero_mat(:,:,j), T, B, x0, xT(:,k), rho, S);
            U_z_all(j) = U_z;
        end
        u_zero_pr = [u_zero_pr; mean(U_z_all)];
    end
end

save([R_dir, 'grad/opt_energy_pr.mat'], 'band_order_pr', 'subj_order_pr', 'u_high_pr', 'u_high2_pr', 'u_high3_pr', 'u_low_pr', 'u_zero_pr')

