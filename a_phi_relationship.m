%% Model Validation - Solve for A, check for linear relationship

% we are making assumptions that the regularization makes this suitable to
% think of as a pathway for information to flow

% We will solve for the C(t) that guves y_dot from y(t), and see if it is
% linearly related to phi(t), the functional connectivity

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

nNode = 102;
nEdges = (nNode^2-nNode)/2;
load([save_dir, 'noise_sg.mat']);

% initialize
high_corr = cell(nSubj, numel(bands));
low_corr = cell(nSubj, numel(bands));

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

% relax B
B(~B) = 1e-5;


%% Loop through data

sens = 'grad';

cnt = 1; % for order
c_all = cell(nSubj,1);

for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    c = cell(numel(bands));
    for k = 1:numel(bands)        
        
        f = bands{k};
        fprintf('\nSubj %03d... band %s...', i, f)

        % get labels
        band_order{cnt} = f;
        subj_order{cnt} = subj;
        cnt = cnt + 1;

        % load states
        load([top_dir, 'states/' subj, '_states.mat'])
        
        % initialize time seris of correlation values
        nTrial = size(states{k},2);
        c_curr = zeros(nNode,nNode,nTrial);
        
        for t = 2:nTrial
            x = states{k}(:,(t-1));
            x_dot = states{k}(:,(t)) - x;

            c_curr(:,:,t) = (x_dot*x')/(x*x');
        end
        c{k} = c_curr;
    end
    c_all{i} = c;
end


for t = 1:nTrial
    figure(1);clf
    imagesc(c_curr(:,:,t))
    title(num2str(t))
    colorbar
    pause(0.01)
end
