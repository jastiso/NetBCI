%% Optimal control to attention pattern

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
rho = [.01, .5];
T = [.2, .05];
x0 = zeros(nNode,1);
x0_cent = zeros(nNode,1) + 1;
x0_321 = zeros(nNode,1) + 2;
S = zeros(nNode,nNode,numel(bands));

% tolerance for error
tol = 1e-3;

% initialize
% first alternate parameter
u_high3_1 = [];
u_low_1 = [];
band_order = {};
subj_order = {};
slope = [];
% second alternate parameter
u_high3_2 = [];
u_low_2 = [];
error = [];


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
S(:,:,1) = eye(nNode); 
xT_attend(:,1) = -(rp + lp);

% beta - contralateral suppression, and ipsilateral activation
S(:,:,2) = eye(nNode); 
xT_attend_mag(:,2) = (-vert);

% gamma - contralateral activation
S(:,:,3) = eye(nNode);
xT_attend(:,3) = (-B - B_control) + lf + rf + lo + ro;

% relax B
B(~B) = 1e-5;

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
        
        %3rd
        high3_mat = get_sg_matrix(nNode, subset(bSG3,:));
        high3_mat = high3_mat./eigs(high3_mat,1) - eye(nNode).*1.001;
        %lowest
        low_mat = get_sg_matrix(nNode, subset(nzSG,:));
        low_mat = low_mat./eigs(low_mat,1) - eye(nNode).*1.001;
        

        % optim_fun(A, T, B, x0, xf, rho, S, error_tol)
        % 3rd
        [~, U_h3, ~] = get_opt_energy(high3_mat, T(1), diag(B), x0, xT_attend(:,k), rho(1), S(:,:,k), tol);
        u_high3_1 = [u_high3_1; U_h3];
        [~, U_h3, err] = get_opt_energy(high3_mat, T(2), diag(B), x0, xT_attend(:,k), rho(2), S(:,:,k), tol);
        u_high3_2 = [u_high3_2; U_h3];

        % low
        [~, U_l, ~] = get_opt_energy(low_mat, T(1), diag(B), x0, xT_attend(:,k), rho(1), S(:,:,k), tol);
        u_low_1 = [u_low_1; U_l];
        [~, U_l, ~] = get_opt_energy(low_mat, T(2), diag(B), x0, xT_attend(:,k), rho(2), S(:,:,k), tol);
        u_low_2 = [u_low_2; U_l];

    end
end

save([R_dir, 'grad/opt_energy_params.mat'], 'band_order', 'subj_order', 'slope', 'u_high3_1', 'u_low_1', 'u_high3_2', 'u_low_2')
save([save_dir, 'oc_error_params'], 'error');


