%% Optimal control to motor imagery and attention pattern

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

% load behavior
load([top_dir, 'Behavior/stats'])

nNode = 102;
nEdges = (nNode^2-nNode)/2;
load([save_dir, 'noise_sg.mat']);

% control parameters
rho = .1;
T = 0.1;
x0 = zeros(nNode,1);
x02 = zeros(nNode,1) + 1;
S = zeros(nNode,nNode,numel(bands));

% tolerance for error
tol = 1e-3;

% initialize
u_high = [];
u_high2 = [];
u_high3 = [];
u_low = [];
u_zero = [];
band_order = {};
subj_order = {};
slope = [];
u_high_a = [];
u_high2_a = [];
u_high3_a = [];
u_low_a = [];
u_high_a2 = [];
u_high2_a2 = [];
u_high3_a2 = [];
u_low_a2 = [];
error = [];
error_a = [];
error_c = [];
error_pr = [];

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
load([top_dir, 'montages/Left_parietal_idx.mat'])
lp = idx;
load([top_dir, 'montages/Right_parietal_idx.mat'])
rp = idx;
load([top_dir, 'montages/Vertex_idx.mat'])
vert = idx;
load([top_dir, 'montages/Left_frontal_idx.mat'])
lf = idx;
load([top_dir, 'montages/Right_frontal_idx.mat'])
rf = idx;
load([top_dir, 'montages/Left_occipital_idx.mat'])
lo = idx;
load([top_dir, 'montages/Right_occipital_idx.mat'])
ro = idx;

% alpha, suppression in left motor
xT(:,1) = -B;
S(:,:,1) = eye(nNode); %diag((xT(:,1) ~= 0));
xT_attend(:,1) = -(rp + lp); % suppression in parietal
xT_attend2(:,1) = -(rp + lp).*2 + 1;

% beta - contralateral suppression, and ipsilateral activation
xT(:,2) = -B + B_control;
S(:,:,2) = eye(nNode); %diag((xT(:,2) ~= 0));
xT_attend(:,2) = -vert; % suppression in midline
xT_attend2(:,2) = (-vert).*2 + 1;

% gamma - contralateral activation
xT(:,3) = B;
S(:,:,3) = eye(nNode); %diag((xT(:,2) ~= 0));
xT_attend(:,3) = (-B - B_control) + lf + rf + lo + ro; 
xT_attend2(:,3) = ((-B - B_control) + lf + rf + lo + ro).*2 + 1; 

% relax B
B(~B) = 1e-5;
B_control(~B_control) = 1e-5;

%% Loop through data

sens = 'grad';

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
        slope(cnt) = betas(i);
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
            zero_mat(:,:,j) = get_sg_matrix(nNode, subset(nbSG(j),:));
            zero_mat(:,:,j) = zero_mat(:,:,j)./eigs(zero_mat(:,:,j),1) - eye(nNode).*1.001;
        end

        % optim_fun(A, T, B, x0, xf, rho, S)
        [x_h, U_h, err] = get_opt_energy(high_mat, T, diag(B), x0, xT(:,k), rho, S(:,:,k));
        u_high = [u_high; U_h];
        [x_ha, U_ha, erra] = get_opt_energy(high_mat, T, diag(B), x0, xT_attend(:,k), rho, S(:,:,k));
        u_high_a = [u_high_a; U_ha];
        error = [error, err];
        error_a = [error_a, erra];
        [~, U_ha2, ~] = get_opt_energy(high_mat, T, diag(B), x02, -xT_attend2(:,k), rho, S(:,:,k));
        u_high_a2 = [u_high_a2; U_ha2];
        
        [x_h2, U_h2, err] = get_opt_energy(high2_mat, T, diag(B), x0, xT(:,k), rho, S(:,:,k));
        u_high2 = [u_high2; U_h2];
        [x_h2a, U_h2a, erra] = get_opt_energy(high2_mat, T, diag(B), x0, xT_attend(:,k), rho, S(:,:,k));
        u_high2_a = [u_high2_a; U_h2a];
        error = [error, err];
        error_a = [error_a, erra];
        [~, U_h2a2, ~] = get_opt_energy(high2_mat, T, diag(B), x02, -xT_attend2(:,k), rho, S(:,:,k));
        u_high2_a2 = [u_high2_a2; U_h2a2];
    
        [x_h3, U_h3, err] = get_opt_energy(high3_mat, T, diag(B), x0, xT(:,k), rho, S(:,:,k));
        u_high3 = [u_high3; U_h3];
        [x_h3a, U_h3a, erra] = get_opt_energy(high3_mat, T, diag(B), x0, xT_attend(:,k), rho, S(:,:,k));
        u_high3_a = [u_high3_a; U_h3a];
        error = [error, err];
        error_a = [error_a, erra];
        [~, U_h3a2, ~] = get_opt_energy(high3_mat, T, diag(B), x02, -xT_attend2(:,k), rho, S(:,:,k));
        u_high3_a2 = [u_high3_a2; U_h3a2];
        
        [x_l, U_l, err] = get_opt_energy(low_mat, T, diag(B), x0, xT(:,k), rho, S(:,:,k));
        u_low = [u_low; U_l];
        [x_la, U_l_a, erra] = get_opt_energy(low_mat, T, diag(B), x0, xT_attend(:,k), rho, S(:,:,k));
        u_low_a = [u_low_a; U_l_a];
        error = [error, err];
        error_a = [error_a, erra];
        [~, U_la2, ~] = get_opt_energy(low_mat, T, diag(B), x02, -xT_attend2(:,k), rho, S(:,:,k));
        u_low_a2 = [u_low_a2; U_la2];
        
        U_z_all = zeros(nZero,1);
        for j = 1:nZero
            [x_z, U_z, err] = get_opt_energy(zero_mat(:,:,j), T, diag(B), x0, xT(:,k), rho, S(:,:,k));
            U_z_all(j) = U_z;
        end
        u_zero = [u_zero; mean(U_z_all)];
        error = [error, err];
        
        
        
        
        % check that you are reaching your target state
        state_idx = logical(diag(S(:,:,k)));
        if norm((x_h(state_idx) - xT(state_idx,k))) > tol
            warning ('Your high error is too large and you are not reaching your target state. Check that your error is within your desired tolerance. If not, try adding more entries to B, fewer to S, or using smaller matrices')
        end
        if norm((x_h2(state_idx) - xT(state_idx,k))) > tol
            warning ('Your high2 error is too large and you are not reaching your target state. Check that your error is within your desired tolerance. If not, try adding more entries to B, fewer to S, or using smaller matrices')
        end
        if norm((x_h3(state_idx) - xT(state_idx,k))) > tol
            warning ('Your high3 error is too large and you are not reaching your target state. Check that your error is within your desired tolerance. If not, try adding more entries to B, fewer to S, or using smaller matrices')
        end
        if norm((x_l(state_idx) - xT(state_idx,k))) > tol
            warning ('Your low error is too large and you are not reaching your target state. Check that your error is within your desired tolerance. If not, try adding more entries to B, fewer to S, or using smaller matrices')
        end
        if norm((x_z(state_idx) - xT(state_idx,k))) > tol
            warning ('Your zero error is too large and you are not reaching your target state. Check that your error is within your desired tolerance. If not, try adding more entries to B, fewer to S, or using smaller matrices')
        end
        
        if norm((x_ha(state_idx) - xT_attend(state_idx,k))) > tol
            warning ('Your attn high error is too large and you are not reaching your target state. Check that your error is within your desired tolerance. If not, try adding more entries to B, fewer to S, or using smaller matrices')
        end
        if norm((x_h2a(state_idx) - xT_attend(state_idx,k))) > tol
            warning ('Your attn high2 error is too large and you are not reaching your target state. Check that your error is within your desired tolerance. If not, try adding more entries to B, fewer to S, or using smaller matrices')
        end
        if norm((x_h3a(state_idx) - xT_attend(state_idx,k))) > tol
            warning ('Your attn high3 error is too large and you are not reaching your target state. Check that your error is within your desired tolerance. If not, try adding more entries to B, fewer to S, or using smaller matrices')
        end
        if norm((x_la(state_idx) - xT_attend(state_idx,k))) > tol
            warning ('Your attn low error is too large and you are not reaching your target state. Check that your error is within your desired tolerance. If not, try adding more entries to B, fewer to S, or using smaller matrices')
        end

    end
end

save([R_dir, 'grad/opt_energy.mat'], 'band_order', 'subj_order', 'slope', 'u_high', 'u_high2', 'u_high3', 'u_low', 'u_zero', ...
    'u_high_a', 'u_high2_a', 'u_high3_a', 'u_low_a', 'u_high_a2', 'u_high2_a2', 'u_high3_a2', 'u_low_a2')
save([save_dir, 'oc_error_mi'], 'error');
save([save_dir, 'oc_error_attn'], 'error_a');


%% Loop through data - other control set

% initialize
u_high_c = [];
u_high2_c = [];
u_high3_c = [];
u_low_c = [];
u_zero_c = [];
band_order_c = {};
subj_order_c = {};
slope = [];

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
        slope(cnt) = betas(i);
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
            zero_mat(:,:,j) = get_sg_matrix(nNode, subset(nbSG(j),:));
            zero_mat(:,:,j) = zero_mat(:,:,j)./eigs(zero_mat(:,:,j),1) - eye(nNode).*1.001;
        end

        % optim_fun(A, T, B, x0, xf, rho, S)
        [~, U_h, err_h] = get_opt_energy(high_mat, T, diag(B_control), x0, xT_attend(:,k), rho, S(:,:,k));
        u_high_c = [u_high_c; U_h];
        [~, U_h2, err_h2] = get_opt_energy(high2_mat, T, diag(B_control), x0, xT_attend(:,k), rho, S(:,:,k));
        u_high2_c = [u_high2_c; U_h2];
        [~, U_h3, err_h3] = get_opt_energy(high3_mat, T, diag(B_control), x0, xT_attend(:,k), rho, S(:,:,k));
        u_high3_c = [u_high3_c; U_h3];
        [~, U_l, err_l] = get_opt_energy(low_mat, T, diag(B_control), x0, xT_attend(:,k), rho, S(:,:,k));
        u_low_c = [u_low_c; U_l];
        U_z_all = zeros(nZero,1);
        for j = 1:nZero
            [~, U_z, err_z] = get_opt_energy(zero_mat(:,:,j), T, diag(B_control), x0, xT_attend(:,k), rho, S(:,:,k));
            U_z_all(j) = U_z;
        end
        u_zero_c = [u_zero_c; mean(U_z_all)];
    end
end

save([R_dir, 'grad/opt_energy_control.mat'], 'slope', 'band_order_c', 'subj_order_c', 'u_high_c', 'u_high2_c', 'u_high3_c', 'u_low_c', 'u_zero_c')



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
            zero_mat(:,:,j) = get_sg_matrix(nNode, subset(nbSG(j),:));
            zero_mat(:,:,j) = zero_mat(:,:,j)./eigs(zero_mat(:,:,j),1) - eye(nNode).*1.001;
        end

        % optim_fun(A, T, B, x0, xf, rho, S)
        [~, U_h, err_h] = get_opt_energy(high_mat, T, diag(B), x0, xT_attend(:,k), rho, S(:,:,k));
        u_high_pr = [u_high_pr; U_h];
        [~, U_h2, err_h2] = get_opt_energy(high2_mat, T, diag(B), x0, xT_attend(:,k), rho, S(:,:,k));
        u_high2_pr = [u_high2_pr; U_h2];
        [~, U_h3, err_h3] = get_opt_energy(high3_mat, T, diag(B), x0, xT_attend(:,k), rho, S(:,:,k));
        u_high3_pr = [u_high3_pr; U_h3];
        [~, U_l, err_l] = get_opt_energy(low_mat, T, diag(B), x0, xT_attend(:,k), rho, S(:,:,k));
        u_low_pr = [u_low_pr; U_l];
        U_z_all = zeros(nZero,1);
        for j = 1:nZero
            [~, U_z, err_z] = get_opt_energy(zero_mat(:,:,j), T, diag(B), x0, xT_attend(:,k), rho, S(:,:,k));
            U_z_all(j) = U_z;
        end
        u_zero_pr = [u_zero_pr; mean(U_z_all)];
    end
end

save([R_dir, 'grad/opt_energy_pr.mat'], 'band_order_pr', 'subj_order_pr', 'slope', 'u_high_pr', 'u_high2_pr', 'u_high3_pr', 'u_low_pr', 'u_zero_pr')

